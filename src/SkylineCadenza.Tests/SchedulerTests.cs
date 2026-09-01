using SkylineCadenza.Core.Scheduling;
using Xunit;

namespace SkylineCadenza.Tests;

public class SchedulerTests
{
    private static Candidate Make(
        string id,
        double mz,
        double rtStart,
        double rtStop,
        string group,
        double quantity = 1e6,
        double[]? fragments = null,
        double? predicted = null,
        CandidateSource source = CandidateSource.Experimental,
        double qValue = 0.001)
    {
        // Synthesise FragmentIon entries from the supplied mz array.
        // Intensities are assigned descending so the input order matches
        // the "top-4 by intensity" derivation, and charge defaults to 1.
        var mzList = fragments ?? new[] { 100.0, 200.0, 300.0, 400.0 };
        var ions = new FragmentIon[mzList.Length];
        for (int i = 0; i < mzList.Length; i++)
            ions[i] = new FragmentIon(mzList[i], mzList.Length - i, 1);
        return new Candidate
        {
            PrecursorId = id,
            StrippedSequence = id,
            ModifiedSequence = id,
            PrecursorCharge = 2,
            PrecursorMz = mz,
            RtStart = rtStart,
            RtStop = rtStop,
            RtApex = (rtStart + rtStop) / 2.0,
            PrecursorQuantity = quantity,
            PredictedPriority = predicted,
            QValue = qValue,
            ProteinQValue = 0.001,
            Proteotypic = 1,
            ProteinGroup = group,
            PeptideType = "unique",
            Fragments = ions,
            Top4Fragments = Candidate.DeriveTopMz(ions, 4),
            Run = "test",
            Source = source,
        };
    }

    [Fact]
    public void Single_precursor_in_an_empty_schedule_always_fits()
    {
        var cands = new[] { Make("A", 500.0, 5.0, 5.4, "G1") };
        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 1 });

        Assert.Single(result.ScheduledIndices);
        Assert.Single(result.Slots);
        Assert.Equal(1, result.ProteinGroupsCovered);
    }

    [Fact]
    public void Two_co_eluting_precursors_with_disjoint_fragments_share_one_slot()
    {
        // Same m/z window, overlapping RT (5.1..5.4 intersection), disjoint
        // top-4 fragments. Must share a slot.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("B", 500.5, 5.1, 5.5, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 1 });

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Single(result.Slots);
        Assert.Equal(2, result.Slots[0].MemberIndices.Count);
    }

    [Fact]
    public void Strict_co_elution_blocks_chain_extension()
    {
        // A and B overlap; B and C overlap; A and C do not. Under the
        // bounding-box bug C would chain in. Under the intersection check
        // (this implementation) C must start a fresh slot.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.6, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("B", 500.5, 5.5, 6.1, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
            Make("C", 500.7, 6.0, 6.6, "G3", fragments: new[] { 175.0, 275.0, 375.0, 475.0 }),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 100 });

        // A and B can share, but C must NOT join that slot (it doesn't co-elute
        // with A). C goes into a separate slot.
        Assert.Equal(3, result.ScheduledIndices.Length);
        Assert.True(result.Slots.Length >= 2);

        // Find C's slot - it cannot contain A.
        int cIdx = Array.IndexOf(result.ScheduledIndices, 2);
        int cSlotId = result.ScheduledSlotIds[cIdx];
        var cSlot = result.Slots.First(s => s.Id == cSlotId);
        int aIdx = Array.IndexOf(result.ScheduledIndices, 0);
        Assert.DoesNotContain(0, cSlot.MemberIndices);
    }

    [Fact]
    public void Fragment_clash_prevents_slot_sharing()
    {
        // Same RT window, similar m/z, but fragments overlap. Must not share.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("B", 500.5, 5.1, 5.5, "G2", fragments: new[] { 100.2, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 100, FragmentTolDa = 0.5 });

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Equal(2, result.Slots.Length);
    }

    [Fact]
    public void Mtm_slot_fit_rejects_precursor_whose_solo_window_would_overhang_slot_edge()
    {
        // Slot-edge rule: every member's PrmIsolationWidthTh quadrupole
        // window must fit fully inside the slot. So every member center
        // must sit at least PrmIsolationWidthTh / 2 from each edge.
        //
        // With IsolationWindowTh = 3.0 and PrmIsolationWidthTh = 0.7 the
        // edge clearance is 0.35 Th. Centers 2.4 Th apart cannot both
        // satisfy that: whichever member's edge is closer would sit
        // 0.30 Th from the slot edge, overhanging by 0.05 Th. The
        // scheduler must split them into two slots.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("B", 502.4, 5.1, 5.5, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            IsolationWindowTh = 3.0,
            PrmIsolationWidthTh = 0.7,
        });

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Equal(2, result.Slots.Length);
    }

    [Fact]
    public void Mtm_slot_fit_merges_precursors_when_both_solo_windows_clear_the_slot_edges()
    {
        // Same edge rule as above. Centers 2.2 Th apart leave 0.40 Th
        // of clearance at the tighter edge - 0.05 Th past the 0.35 Th
        // the solo window needs - so both members' quadrupole windows
        // fit inside the slot and the scheduler should multiplex them
        // into one slot.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("B", 502.2, 5.1, 5.5, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            IsolationWindowTh = 3.0,
            PrmIsolationWidthTh = 0.7,
        });

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Single(result.Slots);
        Assert.Equal(2, result.Slots[0].MemberIndices.Count);
    }

    [Fact]
    public void Prm_mode_forces_one_precursor_per_slot()
    {
        // Same inputs that would multiplex in MTM mode; PRM should keep them
        // in separate slots.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("B", 500.5, 5.1, 5.5, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            Mode = AcquisitionMode.Prm,
            CycleBudget = 100,
        });

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Equal(2, result.Slots.Length);
    }

    [Fact]
    public void Budget_of_one_with_non_co_eluting_precursors_drops_the_later_arrival()
    {
        // Two solo precursors at different RTs both want a slot. Budget=1
        // means only one slot can be active at a time. They don't overlap in
        // RT, so a single slot CAN cover both at different times - but they
        // can't share an MTM window because they don't co-elute. With
        // budget=1, the second precursor can still be scheduled into a fresh
        // slot because at any moment only one slot is active. This test
        // therefore expects BOTH to schedule (different RT windows = no
        // concurrency violation).
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1"),
            Make("B", 600.0, 8.0, 8.4, "G2"),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 1 });

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Equal(2, result.Slots.Length);
    }

    [Fact]
    public void Min_peptides_per_protein_drops_groups_below_threshold()
    {
        // Two protein groups. G1 has 2 peptides that can both be scheduled;
        // G2 has only 1 peptide. With MinPeptidesPerProtein = 2, only G1
        // survives the post-filter.
        var cands = new[]
        {
            Make("G1.a", 500.0, 5.0, 5.4, "G1", quantity: 1e8),
            Make("G1.b", 600.0, 6.0, 6.4, "G1", quantity: 1e7),
            Make("G2.a", 700.0, 7.0, 7.4, "G2", quantity: 1e6),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MinPeptidesPerProtein = 2,
            MaxPeptidesPerProtein = 5,
        });

        // G1's two peptides survive; G2's one peptide is dropped.
        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Equal(1, result.ProteinGroupsCovered);

        // The slot belonging to G2 should also be gone.
        Assert.All(result.Slots, s =>
            Assert.All(s.MemberIndices, idx => Assert.Equal("G1", cands[idx].ProteinGroup)));
    }

    [Fact]
    public void Min_peptides_default_of_one_keeps_singletons()
    {
        var cands = new[]
        {
            Make("G1.a", 500.0, 5.0, 5.4, "G1"),
            Make("G2.a", 700.0, 7.0, 7.4, "G2"),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 100 });
        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Equal(2, result.ProteinGroupsCovered);
    }

    [Fact]
    public void MaximizeProteins_prefers_joinable_peptide_to_free_slot_for_another_protein()
    {
        // Budget = 2, four proteins all co-eluting in [5.0, 5.4].
        //
        // G1: A1 mz 500.0                             - opens slot s1.
        // G2: B1 mz 600.0  (best score, would open new slot)
        //     B2 mz 500.5  (lower score, joins s1: a 0.5 Th center
        //                   gap easily clears the slot-edge rule for
        //                   the default 3.0 / 0.7 settings, and the
        //                   fragments don't clash with A1).
        // G3: C  mz 700.0  (needs a new slot).
        // G4: D  mz 800.0  (needs a new slot).
        //
        // Balanced: cover picks G2's best (B1) → opens slot s2, slot
        // count at this RT = 2. G3 can't open a third slot → dropped.
        // G4 dropped. Total 2 proteins.
        //
        // MaximizeProteins: cover sees B2 is joinable, takes it →
        // s1 carries G1+G2 (count stays 1). G3 opens s2 (count 2).
        // G4 still blocked. Total 3 proteins.
        var cands = new[]
        {
            Make("G1.A1", 500.0, 5.0, 5.4, "G1", quantity: 5e6,
                fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("G2.B1", 600.0, 5.0, 5.4, "G2", quantity: 3e6,
                fragments: new[] { 110.0, 210.0, 310.0, 410.0 }),
            Make("G2.B2", 500.5, 5.0, 5.4, "G2", quantity: 1e6,
                fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
            Make("G3.C",  700.0, 5.0, 5.4, "G3", quantity: 3e6,
                fragments: new[] { 120.0, 220.0, 320.0, 420.0 }),
            Make("G4.D",  800.0, 5.0, 5.4, "G4", quantity: 2e6,
                fragments: new[] { 130.0, 230.0, 330.0, 430.0 }),
        };

        var balanced = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 2,
            Objective = CoverageObjective.Balanced,
        });
        var maxProteins = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 2,
            Objective = CoverageObjective.MaximizeProteins,
        });

        Assert.Equal(2, balanced.ProteinGroupsCovered);
        Assert.Equal(3, maxProteins.ProteinGroupsCovered);
    }

    [Fact]
    public void MaximizePeptides_uncaps_load_up_but_keeps_round_robin_fairness()
    {
        // Two proteins, four peptides each, well-separated RTs so every
        // peptide opens its own slot (no MTM merging here). With Max = 2:
        //
        //   Balanced: 2 peptides per protein → 4 total.
        //   MaximizePeptides: load-up keeps going past Max →
        //                     all 8 peptides scheduled.
        //
        // Round-robin order is preserved by construction: the load-up
        // loop iterates groups in order and adds at most one peptide per
        // group per pass.
        var cands = new[]
        {
            Make("G1.a", 500.0, 5.0, 5.4, "G1", quantity: 4e6),
            Make("G1.b", 510.0, 6.0, 6.4, "G1", quantity: 3e6),
            Make("G1.c", 520.0, 7.0, 7.4, "G1", quantity: 2e6),
            Make("G1.d", 530.0, 8.0, 8.4, "G1", quantity: 1e6),
            Make("G2.a", 600.0, 5.5, 5.9, "G2", quantity: 4e6),
            Make("G2.b", 610.0, 6.5, 6.9, "G2", quantity: 3e6),
            Make("G2.c", 620.0, 7.5, 7.9, "G2", quantity: 2e6),
            Make("G2.d", 630.0, 8.5, 8.9, "G2", quantity: 1e6),
        };

        var balanced = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 2,
            Objective = CoverageObjective.Balanced,
        });
        var maxPeptides = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 2,
            Objective = CoverageObjective.MaximizePeptides,
        });

        Assert.Equal(4, balanced.ScheduledIndices.Length);
        Assert.Equal(8, maxPeptides.ScheduledIndices.Length);

        // Fairness check: at every prefix of the load-up sequence the
        // per-group counts differ by at most 1. We can verify the final
        // state symmetrically: both groups end with 4 peptides each.
        int g1 = maxPeptides.ScheduledIndices.Count(i => cands[i].ProteinGroup == "G1");
        int g2 = maxPeptides.ScheduledIndices.Count(i => cands[i].ProteinGroup == "G2");
        Assert.Equal(4, g1);
        Assert.Equal(4, g2);
    }

    [Fact]
    public void Exclusive_target_list_matches_bare_accession_against_UniProt_style_group()
    {
        // Skyline-document ingest carries the full Skyline protein name
        // (e.g. "sp|P55011|S12A2_HUMAN") on each candidate's
        // ProteinGroup, while ProteinListParser extracts bare accessions
        // ("P55011") from FASTA target lists. Exclusive mode used to
        // exact-string Contains, dropping every UniProt-style group
        // silently; the filter now normalises one side so the
        // intersection actually happens.
        var cands = new[]
        {
            Make("PEP1", 500.0, 5.0, 5.4, "sp|P55011|S12A2_HUMAN"),
            Make("PEP2", 600.0, 6.0, 6.4, "sp|P55011|S12A2_HUMAN"),
            Make("PEP3", 700.0, 7.0, 7.4, "sp|Q9NZW4|DSPP_HUMAN"),
            Make("PEP4", 800.0, 8.0, 8.4, "sp|P12345|FOOBAR_HUMAN"),
        };
        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            TargetProteins = new HashSet<string> { "P55011", "Q9NZW4" },
            TargetMode = TargetListMode.Exclusive,
        });

        // P12345 should be excluded; the other two should be present.
        var groupsCovered = result.ScheduledIndices
            .Select(i => cands[i].ProteinGroup).Distinct().ToHashSet();
        Assert.Equal(2, groupsCovered.Count);
        Assert.Contains("sp|P55011|S12A2_HUMAN", groupsCovered);
        Assert.Contains("sp|Q9NZW4|DSPP_HUMAN", groupsCovered);
        Assert.DoesNotContain("sp|P12345|FOOBAR_HUMAN", groupsCovered);
    }

    [Fact]
    public void Exclusive_target_list_still_works_when_group_is_bare_accession()
    {
        // Regression check: don't break the existing path where the
        // candidate's ProteinGroup is already a bare accession (DIA-NN
        // path, Carafe path).
        var cands = new[]
        {
            Make("PEP1", 500.0, 5.0, 5.4, "P55011"),
            Make("PEP2", 600.0, 6.0, 6.4, "P12345"),
        };
        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            TargetProteins = new HashSet<string> { "P55011" },
            TargetMode = TargetListMode.Exclusive,
        });
        var groupsCovered = result.ScheduledIndices
            .Select(i => cands[i].ProteinGroup).Distinct().ToHashSet();
        Assert.Single(groupsCovered);
        Assert.Contains("P55011", groupsCovered);
    }

    [Fact]
    public void Load_up_fills_to_Max_under_both_Balanced_and_MaximizeProteins()
    {
        // One protein, three candidates that don't overlap in RT. Cover
        // pass schedules the highest-score one; the load-up loop then
        // adds the remaining two up to MaxPeptidesPerProtein (default
        // 5) as long as the cycle budget allows.
        //
        // Both Balanced and MaximizeProteins cap load-up at Max, so
        // both should end up scheduling all three peptides here. The
        // difference between these objectives is in the cover-pass
        // strategy (MaximizeProteins prefers peptides that join an
        // existing slot), not in the load-up cap.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1", quantity: 1e8),
            Make("B", 600.0, 6.0, 6.4, "G1", quantity: 1e7),
            Make("C", 700.0, 7.0, 7.4, "G1", quantity: 1e6),
        };

        var balanced = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            Objective = CoverageObjective.Balanced,
        });
        var maxProteins = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            Objective = CoverageObjective.MaximizeProteins,
        });

        Assert.Equal(3, balanced.ScheduledIndices.Length);
        Assert.Equal(3, maxProteins.ScheduledIndices.Length);
    }

    [Fact]
    public void Load_up_distributes_peptides_breadth_first_across_proteins()
    {
        // Five proteins, three candidates each, all at distinct
        // non-overlapping RTs. Generous cycle budget so every peptide
        // could fit. After cover pass each protein has 1; the
        // breadth-first round-robin load-up should bring every
        // protein up to Max = 3 in lockstep, not load any single
        // protein to Max before others get their 2nd peptide.
        var cands = new List<Candidate>();
        string[] proteins = { "G1", "G2", "G3", "G4", "G5" };
        double rt = 5.0;
        foreach (var g in proteins)
        {
            for (int j = 0; j < 3; j++)
            {
                cands.Add(Make(
                    id: $"{g}_p{j}",
                    mz: 400.0 + rt * 10,
                    rtStart: rt,
                    rtStop: rt + 0.3,
                    group: g,
                    quantity: 1e9 - rt));
                rt += 0.8;
            }
        }

        var result = Scheduler.Run(cands.ToArray(), new SchedulingParameters
        {
            CycleBudget = 100,
            Objective = CoverageObjective.MaximizeProteins,
            MaxPeptidesPerProtein = 3,
        });

        var perGroup = result.ScheduledIndices
            .GroupBy(i => cands[i].ProteinGroup)
            .ToDictionary(g => g.Key, g => g.Count());

        Assert.Equal(proteins.Length, perGroup.Count);
        foreach (var g in proteins)
            Assert.Equal(3, perGroup[g]);
        Assert.Equal(15, result.ScheduledIndices.Length);
    }

    [Fact]
    public void MaximizeProteins_load_up_honours_MaxPeptidesPerProtein()
    {
        // One protein, five candidates that don't overlap in RT and
        // would all fit individually inside the cycle budget. The cover
        // pass takes the best-scoring one; load-up should fill the rest
        // up to Max = 3, not all 5.
        var cands = new[]
        {
            Make("A", 500.0, 5.0, 5.4, "G1", quantity: 5e8),
            Make("B", 600.0, 6.0, 6.4, "G1", quantity: 4e8),
            Make("C", 700.0, 7.0, 7.4, "G1", quantity: 3e8),
            Make("D", 800.0, 8.0, 8.4, "G1", quantity: 2e8),
            Make("E", 900.0, 9.0, 9.4, "G1", quantity: 1e8),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            Objective = CoverageObjective.MaximizeProteins,
            MaxPeptidesPerProtein = 3,
        });

        Assert.Equal(3, result.ScheduledIndices.Length);
    }

    [Fact]
    public void Scheduler_PredictedRanking_OrdersPeptidesWithinProteinByScoreDescNullsLast()
    {
        // Three peptides in one protein at distinct, non-overlapping RTs. With
        // Max = 1 only the top-priority peptide is covered. PredictedRanking
        // must pick the highest predicted score; a null-scored peptide (no
        // overlay match) sorts last and must never win.
        var cands = new[]
        {
            Make("low", 500.0, 1.0, 1.4, "G1", predicted: 0.2),
            Make("high", 600.0, 5.0, 5.4, "G1", predicted: 0.9),
            Make("none", 700.0, 9.0, 9.4, "G1", predicted: null),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MinPeptidesPerProtein = 1,
            MaxPeptidesPerProtein = 1,
            PeptideRanking = PeptidePriority.PredictedRanking,
            Objective = CoverageObjective.MaximizeProteins,
        });

        var scheduledIds = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToList();
        Assert.Single(scheduledIds);
        Assert.Equal("high", scheduledIds[0]);
    }

    [Fact]
    public void Scheduler_ExperimentalWinsScarceBudgetOverPredicted()
    {
        // One concurrent slot available at this RT. A confirmed target and a
        // predicted target both want it but their fragments clash, so they
        // cannot multiplex. The experimental tier runs first and takes the
        // only slot; the predicted target is dropped rather than displacing it.
        var cands = new[]
        {
            Make("exp", 500.0, 5.0, 5.4, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 },
                source: CandidateSource.Experimental),
            Make("pred", 500.2, 5.0, 5.4, "G2", fragments: new[] { 100.1, 200.1, 300.1, 400.1 },
                source: CandidateSource.Predicted),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 1, FragmentTolDa = 0.5 });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToList();
        Assert.Single(scheduled);
        Assert.Equal("exp", scheduled[0]);
    }

    [Fact]
    public void Scheduler_PredictedBackfillsIntoLeftoverSlotCapacity()
    {
        // A confirmed target opens a slot; a co-eluting predicted target with
        // disjoint fragments backfills the same slot (leftover multiplexing
        // capacity) instead of being ignored.
        var cands = new[]
        {
            Make("exp", 500.0, 5.0, 5.4, "G1", fragments: new[] { 100.0, 200.0, 300.0, 400.0 },
                source: CandidateSource.Experimental),
            Make("pred", 500.5, 5.1, 5.5, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 },
                source: CandidateSource.Predicted),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 1 });

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Single(result.Slots);
        Assert.Equal(2, result.Slots[0].MemberIndices.Count);
        Assert.Equal(2, result.ProteinGroupsCovered);
    }

    [Fact]
    public void Scheduler_PredictedRespectsPerProteinCapAlreadySpentByExperimental()
    {
        // Shared per-protein cap across tiers. G1 spends its Max=2 on two
        // confirmed peptides, so its predicted peptide must NOT be added. G2
        // has only one confirmed peptide, so a predicted peptide DOES backfill
        // its remaining depth. All at distinct RTs so budget is never binding.
        var cands = new[]
        {
            Make("G1.exp1", 500.0, 1.0, 1.4, "G1", quantity: 9e6, source: CandidateSource.Experimental),
            Make("G1.exp2", 510.0, 2.0, 2.4, "G1", quantity: 8e6, source: CandidateSource.Experimental),
            Make("G1.pred", 520.0, 3.0, 3.4, "G1", quantity: 7e6, source: CandidateSource.Predicted),
            Make("G2.exp1", 600.0, 5.0, 5.4, "G2", quantity: 9e6, source: CandidateSource.Experimental),
            Make("G2.pred", 610.0, 6.0, 6.4, "G2", quantity: 7e6, source: CandidateSource.Predicted),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 2,
            Objective = CoverageObjective.MaximizeProteins,
        });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToHashSet();
        Assert.Contains("G1.exp1", scheduled);
        Assert.Contains("G1.exp2", scheduled);
        Assert.DoesNotContain("G1.pred", scheduled); // G1 already at cap via experimental
        Assert.Contains("G2.exp1", scheduled);
        Assert.Contains("G2.pred", scheduled);       // G2 had room for predicted depth
    }

    [Fact]
    public void Scheduler_PredictedCoversProteinNotSeenExperimentally()
    {
        // A protein present only in the predicted tier is covered in the
        // second pass once the experimental tier is done.
        var cands = new[]
        {
            Make("G1.exp", 500.0, 5.0, 5.4, "G1", source: CandidateSource.Experimental),
            Make("G2.pred", 800.0, 8.0, 8.4, "G2", source: CandidateSource.Predicted),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 100 });

        Assert.Equal(2, result.ProteinGroupsCovered);
        var scheduled = result.ScheduledIndices.Select(i => cands[i]).ToList();
        Assert.Contains(scheduled, c => c.ProteinGroup == "G2" && c.Source == CandidateSource.Predicted);
    }

    [Fact]
    public void Scheduler_ExperimentalCoreUnchangedWhenBackfillReordersProteins()
    {
        // The experimental tier must schedule identically whether or not a
        // predicted backfill tier is present. Here G1's confirmed peptide
        // outweighs G2's confirmed peptide, so with only one concurrent slot G1
        // wins coverage. G2 also carries a high-intensity PREDICTED peptide: if
        // predicted candidates were allowed to inflate G2's summed intensity, the
        // protein order would flip and G2 would steal the slot - overwriting the
        // empirical schedule. The lock keeps the experimental core fixed.
        var expOnly = new[]
        {
            Make("G1.exp", 500.0, 5.0, 5.4, "G1", quantity: 9e6),
            Make("G2.exp", 700.0, 5.0, 5.4, "G2", quantity: 6e6),
        };
        var withBackfill = new[]
        {
            Make("G1.exp", 500.0, 5.0, 5.4, "G1", quantity: 9e6),
            Make("G2.exp", 700.0, 5.0, 5.4, "G2", quantity: 6e6),
            Make("G2.pred", 900.0, 5.0, 5.4, "G2", quantity: 20e6,
                source: CandidateSource.Predicted),
        };

        var p = new SchedulingParameters { CycleBudget = 1 };
        var a = Scheduler.Run(expOnly, p);
        var b = Scheduler.Run(withBackfill, p);

        string[] ExpCore(Candidate[] cands, ScheduleResult r) =>
            r.ScheduledIndices
                .Where(i => cands[i].Source == CandidateSource.Experimental)
                .Select(i => cands[i].PrecursorId)
                .OrderBy(s => s)
                .ToArray();

        var coreA = ExpCore(expOnly, a);
        var coreB = ExpCore(withBackfill, b);
        Assert.Equal(new[] { "G1.exp" }, coreA);
        Assert.Equal(coreA, coreB); // backfill did not overwrite the empirical core
    }

    [Fact]
    public void Scheduler_ConfirmedTierIgnoresPredictedRankingMetric()
    {
        // The app flips PeptideRanking to PredictedRanking whenever a Carafe
        // detectability / Bromo overlay is loaded with backfill. Confirmed
        // peptides carry no PredictedPriority, so ranking THEM by it would
        // collapse them to a PrecursorId tie and pick a different top-N under the
        // per-protein cap. The confirmed tier must keep ranking by
        // PrecursorQuantity so its schedule matches a backfill-off run.
        //
        // exp3 has the lexicographically-smallest PrecursorId but the lowest
        // quantity: under a raw PredictedRanking sort (all confirmed tie at
        // null, tie-broken by PrecursorId) exp3 would win the cap over exp1;
        // under quantity it must not.
        var expOnly = new[]
        {
            Make("Zexp1", 500.0, 1.0, 1.4, "G1", quantity: 9e6),
            Make("Yexp2", 510.0, 2.0, 2.4, "G1", quantity: 8e6),
            Make("Aexp3", 520.0, 3.0, 3.4, "G1", quantity: 1e6),
        };
        var combined = new[]
        {
            Make("Zexp1", 500.0, 1.0, 1.4, "G1", quantity: 9e6),
            Make("Yexp2", 510.0, 2.0, 2.4, "G1", quantity: 8e6),
            Make("Aexp3", 520.0, 3.0, 3.4, "G1", quantity: 1e6),
            Make("G1pred", 530.0, 4.0, 4.4, "G1", quantity: 5e6, predicted: 0.99,
                source: CandidateSource.Predicted),
        };

        var off = new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 2,
            PeptideRanking = PeptidePriority.PrecursorQuantity, // backfill-off default
        };
        var on = off with { PeptideRanking = PeptidePriority.PredictedRanking }; // app's backfill flip

        var rOff = Scheduler.Run(expOnly, off);
        var rOn = Scheduler.Run(combined, on);

        var coreOff = rOff.ScheduledIndices.Select(i => expOnly[i].PrecursorId)
            .OrderBy(s => s).ToArray();
        var coreOn = rOn.ScheduledIndices
            .Where(i => combined[i].Source == CandidateSource.Experimental)
            .Select(i => combined[i].PrecursorId).OrderBy(s => s).ToArray();

        Assert.Equal(new[] { "Yexp2", "Zexp1" }, coreOff); // top-2 by quantity, not exp3
        Assert.Equal(coreOff, coreOn); // backfill flip to PredictedRanking did not move the core
    }

    [Fact]
    public void Scheduler_ConfirmedTierKeepsUserRankingIndependentOfPredictedTier()
    {
        // Enabling backfill must not change the confirmed peptides' ranking. The
        // user ranks confirmed peptides by QValue, which disagrees with intensity:
        // exp_lowq has the best (smallest) q-value but low intensity. With Max=1
        // the confirmed tier keeps exp_lowq (the QValue winner), NOT exp_hiq (the
        // intensity winner) - proving the confirmed tier ranks by the user's
        // metric, not by PredictedRanking or a PrecursorQuantity fallback. The
        // predicted tier meanwhile backfills, ranked by its own detectability.
        var cands = new[]
        {
            Make("G1.exp_lowq", 500.0, 1.0, 1.4, "G1", quantity: 1e6, qValue: 0.0001),
            Make("G1.exp_hiq",  510.0, 2.0, 2.4, "G1", quantity: 9e6, qValue: 0.05),
            Make("G2.pred", 800.0, 8.0, 8.4, "G2", predicted: 0.99,
                source: CandidateSource.Predicted),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 1,
            PeptideRanking = PeptidePriority.QValue,
        });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToHashSet();
        Assert.Contains("G1.exp_lowq", scheduled);       // QValue winner kept
        Assert.DoesNotContain("G1.exp_hiq", scheduled);  // not the intensity winner
        Assert.Contains("G2.pred", scheduled);           // predicted still backfills
    }

    [Fact]
    public void Scheduler_PredictedCannotShrinkOutExperimentalMember()
    {
        // Disjoint packing: a predicted precursor sits within a quad width of a
        // multiplexed experimental slot and cannot open its own window (its band
        // overlaps the slot's) or merge (its fragments clash). The only way in
        // would be to shrink the slot and evict the co-isolated experimental
        // member - which the locked experimental tier forbids. The predicted
        // target is dropped; both experimental members stay put.
        var cands = new[]
        {
            Make("exp1", 500.0, 5.0, 5.4, "G1", quantity: 9e6,
                fragments: new[] { 100.0, 200.0, 300.0, 400.0 }),
            Make("exp2", 501.0, 5.0, 5.4, "G1", quantity: 8e6,
                fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
            Make("pred", 500.2, 5.0, 5.4, "G2", quantity: 20e6,
                fragments: new[] { 100.0, 200.0, 300.0, 400.0 },
                source: CandidateSource.Predicted),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 5,
            FragmentTolDa = 0.5,
            Packing = WindowPacking.Disjoint,
        });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToHashSet();
        Assert.Contains("exp1", scheduled);
        Assert.Contains("exp2", scheduled);
        Assert.DoesNotContain("pred", scheduled);
    }

    [Fact]
    public void Scheduler_PeptideRankCap_LimitsGroupToTopRankedByQuantity()
    {
        // Four peptides in one protein, all at distinct RTs and m/z so budget
        // and geometry never bind. With a rank cap of 2 only the two
        // highest-quantity ones are ELIGIBLE, even though Max would allow 5.
        var cands = new[]
        {
            Make("p1", 500.0, 1.0, 1.4, "G1", quantity: 9e6),
            Make("p2", 510.0, 2.0, 2.4, "G1", quantity: 8e6),
            Make("p3", 520.0, 3.0, 3.4, "G1", quantity: 7e6),
            Make("p4", 530.0, 4.0, 4.4, "G1", quantity: 6e6),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 5,
            PeptideRankCap = 2,
            PeptideRanking = PeptidePriority.PrecursorQuantity,
        });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToHashSet();
        Assert.Equal(new HashSet<string> { "p1", "p2" }, scheduled);
    }

    [Fact]
    public void Scheduler_PeptideRankCap_FollowsQValueRanking()
    {
        // Same pool, but ranked by q-value ascending - which is the exact
        // REVERSE of the quantity order. The cap must keep the top two by
        // q-value, proving it tracks whichever metric PeptideRanking selects.
        var cands = new[]
        {
            Make("p1", 500.0, 1.0, 1.4, "G1", quantity: 9e6, qValue: 0.009),
            Make("p2", 510.0, 2.0, 2.4, "G1", quantity: 8e6, qValue: 0.007),
            Make("p3", 520.0, 3.0, 3.4, "G1", quantity: 7e6, qValue: 0.002),
            Make("p4", 530.0, 4.0, 4.4, "G1", quantity: 6e6, qValue: 0.001),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 5,
            PeptideRankCap = 2,
            PeptideRanking = PeptidePriority.QValue,
        });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToHashSet();
        Assert.Equal(new HashSet<string> { "p3", "p4" }, scheduled);
    }

    [Fact]
    public void Scheduler_PeptideRankCap_FollowsPredictedRanking()
    {
        // Predicted ranking (Carafe detectability / Bromo score), descending.
        // Quantity order is again the reverse, so only the predicted score can
        // produce this result.
        var cands = new[]
        {
            Make("p1", 500.0, 1.0, 1.4, "G1", quantity: 9e6, predicted: 0.10),
            Make("p2", 510.0, 2.0, 2.4, "G1", quantity: 8e6, predicted: 0.40),
            Make("p3", 520.0, 3.0, 3.4, "G1", quantity: 7e6, predicted: 0.95),
            Make("p4", 530.0, 4.0, 4.4, "G1", quantity: 6e6, predicted: 0.99),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 5,
            PeptideRankCap = 2,
            PeptideRanking = PeptidePriority.PredictedRanking,
        });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToHashSet();
        Assert.Equal(new HashSet<string> { "p3", "p4" }, scheduled);
    }

    [Fact]
    public void Scheduler_PeptideRankCap_AppliesPerSourceTier()
    {
        // Backfill pool. A cap of 1 must mean "top 1 confirmed AND top 1
        // predicted", not "top 1 of the merged pool" - the tiers rank by
        // different metrics and the confirmed tier must not spend the
        // predicted tier's cap (or vice versa).
        var cands = new[]
        {
            Make("exp1", 500.0, 1.0, 1.4, "G1", quantity: 9e6),
            Make("exp2", 510.0, 2.0, 2.4, "G1", quantity: 8e6),
            Make("pred1", 520.0, 3.0, 3.4, "G1", quantity: 1e6, predicted: 0.99,
                source: CandidateSource.Predicted),
            Make("pred2", 530.0, 4.0, 4.4, "G1", quantity: 1e6, predicted: 0.80,
                source: CandidateSource.Predicted),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 5,
            PeptideRankCap = 1,
            PeptideRanking = PeptidePriority.PrecursorQuantity,
        });

        var scheduled = result.ScheduledIndices.Select(i => cands[i].PrecursorId).ToHashSet();
        Assert.Equal(new HashSet<string> { "exp1", "pred1" }, scheduled);
    }

    [Fact]
    public void Scheduler_PeptideRankCapZero_SchedulesEveryEligiblePeptide()
    {
        // The default (0 = uncapped) must leave the pool untouched, which is
        // what keeps existing schedules and the golden fixtures stable.
        var cands = new[]
        {
            Make("p1", 500.0, 1.0, 1.4, "G1", quantity: 9e6),
            Make("p2", 510.0, 2.0, 2.4, "G1", quantity: 8e6),
            Make("p3", 520.0, 3.0, 3.4, "G1", quantity: 7e6),
        };

        var result = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MaxPeptidesPerProtein = 5,
            PeptideRankCap = 0,
        });

        Assert.Equal(3, result.ScheduledIndices.Length);
    }
}
