# Scheduler bug: window expansion co-isolates unchecked precursors

Status: **fixed.** Slot-aware packing is the only behaviour; the `WindowPacking`
toggle and its `Legacy` mode were removed. Covered by
[`SchedulerWindowScopingTests.cs`](../src/SkylineCadenza.Tests/SchedulerWindowScopingTests.cs),
which asserts the invariant holds, that only physically unseparable
co-isolations survive, that concurrent bands never overlap, and that PRM is
unaffected.

## Summary

`Scheduler.TrySchedule` validates the wrong set. Its fragment-clash gate
compares a candidate only against the slot's current **members**
(`slot.Fragments`, the union of member top-4 m/z; see
`Scheduler.cs` merge path). But the physical constraint is the isolation
**window**: a quadrupole isolates the entire contiguous band
`[MzMin - PrmIsolationWidthTh/2, MzMax + PrmIsolationWidthTh/2]`, so it
co-isolates every precursor whose m/z falls in that band, member or not.

Because a successful merge widens the band (`MzMin = min(...)`,
`MzMax = max(...)`), reaching a distant precursor silently pulls in everything
between it and the existing members. Nothing re-checks the widened band against
non-member precursors. The result: the scheduler can **explicitly reject a
precursor for clashing and then physically re-include it** by expanding the
window across its m/z. The exclusion is illusory.

## The A / Q / P example

Three co-eluting precursors (same RT window), with the default
`IsolationWindowTh = 3.0` and `PrmIsolationWidthTh = 0.7` (center-span budget
2.3 Th):

| precursor | m/z | fragments (m/z) | clashes with |
|-----------|-----|-----------------|--------------|
| A | 600.00 | 1000, 1100, 1200, 1300 | Q |
| Q | 600.75 | **1000**, **1400**, 1500, 1600 | A and P |
| P | 601.50 | **1400**, 1700, 1800, 1900 | Q |

A and P share no fragment, so they may multiplex; Q shares a fragment with each.

What the scheduler does today:

1. Q is scheduled first (smallest group / covered first) and opens a slot at 600.75.
2. A tries to join Q's slot, clashes on 1000, is **rejected**, and opens its own slot at 600.00.
3. P tries Q's slot, clashes on 1400, is **rejected**; it then joins A's slot
   (span 1.5 Th &le; 2.3, no A/P clash). A's slot window becomes `[600.0, 601.5]`,
   emitted band `[599.65, 601.85]`.

Q (600.75) now sits in the middle of A+P's band. When that slot fires, the quad
co-isolates Q with A, reintroducing exactly the fragment interference (on m/z
1000) that Q was rejected to avoid. The scheduler believes the slot is clean
because it only ever compared A&harr;P.

## Two consequences

- **The exclusion is a no-op.** Rejecting a clashing precursor accomplishes
  nothing if a later expansion crosses its m/z; the instrument does not know it
  "isn't a member".
- **The within-slot guarantee is weaker than advertised.** "No two members
  clash" does not imply "no clash inside the window": non-members (rejected
  candidates, other slots' members, shared/razor peptides, and the un-modelled
  non-target background) live in the same band and are never checked.

The one bound on the blast radius is the center-span cap
(`&le; IsolationWindowTh - PrmIsolationWidthTh`, ~2.3 Th). It stops unbounded
growth but leaves plenty of co-eluting proteome unchecked, and at
whole-proteome saturation that band is dense.

## The invariant we want

No scheduled slot's isolation band may contain a co-eluting precursor that
shares a fragment (within `FragmentTolDa`) with one of that slot's members,
unless that precursor is itself a member and was therefore interference-checked.

## How the fix relates to the redesign

This is the concrete motivation for making the interference check
**window-scoped, not member-scoped** (the "feature 2 / feature 3" rework):

- The clean-fragment / clash check must run over every precursor whose m/z lands
  in `[isoLo, isoHi]` at the slot's RT, drawn from the full candidate pool (and,
  ideally, a non-target background index), not just the member list.
- Reframed, expanding the window is not "add member P" but "commit to
  co-isolating the whole `[MzMin, MzMax]` interval", so the merge test becomes
  "does everything now in that interval still keep its clean fragments". That
  also makes the partition boundary fall exactly where widening would start
  swallowing interferers &mdash; the m/z-sweep window construction described in
  the partition redesign.

The member-scoped check no longer ships: windows are always packed disjointly.

## What slot-aware packing actually does

Two guards, with different lifetimes:

1. **Growth bounds** (`Slot.GrowLo` / `Slot.GrowHi`), tightened by a bounded scan
   of the candidate **pool** (not of the slots built so far) at slot creation and
   on every member add. A window can therefore never expand across a co-eluting
   precursor that clashes with one of its members. Because the pool is static and
   bounds only ever tighten, the result is independent of the order proteins are
   processed in - which is what lets this work inside the existing
   protein-priority greedy.
2. **Neighbour-band probe** at expansion and slot-open time, via the existing
   `slotByMzBin` index: reject when the resulting band would overlap a
   concurrently-firing slot's band. Slots appear as scheduling proceeds, so this
   one is probed rather than cached.

Both are gated to MTM. In PRM a slot holds exactly one precursor and its band is
a fixed `PrmIsolationWidthTh`, so it never grows (the failure is impossible) and
the non-overlap rule could only ever drop targets within a quad width of each
other, which no window scheme can separate.

## Measured cost

On balanced_GPF (4,698 precursors, 3,075 proteins, budget 100), against the
member-scoped behaviour that used to be the default:

| | member-scoped (removed) | slot-aware (current) |
|---|---|---|
| proteins covered | 3,075 | ~2,230 |
| windows | 1,818 | ~1,300 |
| peak concurrency | 100 | ~67 |
| mean isolation width | 1.70 Th | 1.26 Th |
| overlapping concurrent band pairs | 2,144-2,219 | **0** |

Confirmed on the 08/03 HeLa GPF report (70,970 precursors, budget 100): 9,079
targets / 5,525 proteins under the old member-scoped default against 6,476 /
4,731 slot-aware, a 29% target and 14% protein reduction. The signal bought with
it is large: the median target's share of its window's co-isolated ion current
rises from 3.30% to 8.48%, and the fraction of windows where a target holds
under 1% of the fill falls from 12.9% to 6.7%. See `MinTargetIonShare` in
`CLAUDE.md` for the measurement.

Overlap is eliminated outright and windows get materially narrower, at a cost of
roughly 28% of protein coverage. The tell is peak concurrency falling well below
the cycle budget: the schedule stops being budget-limited and becomes
geometry-limited, with no mechanism to spend the leftover capacity. Tightening
`FragmentTolDa` from 0.5 to 0.02 Da recovers only ~25 proteins, so the cost comes
from the blanket neighbour-band rule rather than from over-eager clash barriers.

**Residual co-isolations are expected and irreducible.** At 0.02 Da tolerance,
10 survive; every one has the interferer within 0.7 Th of the member it clashes
with (measured gaps: 0.006-0.145 Th). Two co-eluting precursors closer than the
quadrupole width cannot be given disjoint bands by any window scheme, so the only
remaining lever for those is to drop the lower-priority one.

A **`BoundsOnly`** middle setting - keep the growth bounds, drop the blanket
neighbour-band rule - is the obvious next experiment: it would still fix the
sweep-over-interferer bug while allowing benign (non-clashing) overlaps, and
should recover most of the lost coverage.
