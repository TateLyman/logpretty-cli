# Lysosome / MTORC1 mechanism audit

## Correction to stage 28

Stage 28 presented the bafilomycin result as the strongest phenotype-first hit and recorded proliferation as unknown. Reading the full text changes that materially. The paper's own **figure title** is:

> *"Bafilomycin A1 promotes differentiation, elevates cell death and decreases chondrocyte proliferation in cultured metatarsal bones."*

And the authors state plainly:

> *"the observed growth stimulation was entirely attributed to the promoted chondrocyte hypertrophy without any contribution from cell proliferation or survival."*

So the length gain is real and well measured, but it is produced by bigger terminal cells **while proliferation falls and apoptosis rises**. Under this project's own interpretation rules that is a trade-off, not productive growth. Stage 28's framing was too favourable and is superseded here.

## What was measured (source-derived)

| endpoint | value |
|---|---|
| bone-length gain, 5 d | Baf 880 +/- 28 um; CQ 710 +/- 18 um; p<0.001; n=7 animals (21 bones) Baf, n=13 animals (39 bones) CQ |
| bafilomycin concentration/schedule | 8 nM, continuous in organ culture, measured over 5-6 d |
| chloroquine concentration/schedule | 30 uM, continuous in organ culture |
| IGF1 positive control | 100 ng/ml; Baf and CQ stimulated growth 'to the same extent as' IGF1 |
| proliferation | DECREASED chondrocyte proliferation |
| apoptosis | ELEVATED cell death by TUNEL labelling |
| terminal hypertrophic-cell size | increased, p<0.01, n=5, 25 measurements per bone in mid-hypertrophic zone of the distal growth plate |
| hypertrophic zone | increased hypertrophy on histology |
| proteoglycan / matrix | Safranin O showed NO change in proteoglycan levels |
| p-RPS6 | Baf 10.1 +/- 2.0 fold (n=3, p=0.011); CQ 2.2 +/- 0.16 fold; CQ increase not significant in one analysis (2.0 +/- 0.4 fold, p=0.23, n=2) |
| p-MTOR | 2.0 +/- 1.2 fold, p=0.49 (NOT significant); total MTOR 1.3 +/- 0.3 fold, p=0.49 |
| p-RPS6KB1 (S6K) | in C5.18 chondrocytes: 14 +/- 44% DECREASE, p=0.78 (not significant) |
| EIF4EBP1 | Torin1 dose-dependently decreased RPS6 and EIF4EBP1 phosphorylation in C5.18 chondrocytes |
| MTORC2 readouts | SGK1 phosphorylation unchanged (21+/-35%, p=0.63); p-AKT Ser473 DOWNregulated |
| SQSTM1 / LC3 | SQSTM1 accumulation in chondrocytes (IHC and western); ATG5, SQSTM1, MAP1LC3A assessed by western after 2 d |
| Atg5-cKO response | Baf stimulated growth similarly in ATG5cKO and control metatarsals (n=7 animals/21 bones control, 6 animals/18 bones cKO) -> autophagy-independent |
| Torin1 interaction | Torin1 'attenuated' / 'significantly diminished' the growth-promoting effect of Baf and CQ (n=5-11 animals). ATTENUATION, NOT ABOLITION - necessity was not demonstrated |
| washout / recovery data | NONE. The strings 'washout' and 'recover' do not occur in the full text. No exposure-then-recovery experiment was performed |
| cell-state dependence | Baf DOWN-regulated p-RPS6 in undifferentiated C5.18 cells but UP-regulated it in differentiated chondrocytes |

## What the authors concluded (author interpretation)

| point | statement |
|---|---|
| growth attributed to hypertrophy | 'the observed growth stimulation was entirely attributed to the promoted chondrocyte hypertrophy without any contribution from cell proliferation or survival' |
| authors' own caveat on MTORC1 | 'we think that genetic studies are required to confirm this assumption and to extend it into a physiological setting' |
| authors flag a dose mismatch | 'activation of RPS6 by Baf is 5 times stronger than by CQ whereas the growth-promoting effect is just 24 percent stronger, suggesting that different mechanisms of growth may be involved' |

## What this pipeline infers (not the authors' claim)

| point | inference |
|---|---|
| net phenotype classification | Length gain arises from larger terminal hypertrophic cells while proliferative output falls and apoptosis rises. Under this project's own rules that is a trade-off, not productive growth |
| MTORC1 necessity | Not demonstrated. Torin1 attenuation is consistent with MTORC1 contributing, but Torin1 also inhibits MTORC1 broadly and suppresses growth on its own, so partial attenuation cannot separate necessity from additive suppression |
| durability | Unknown. With no washout arm and a 5-6 d culture, nothing in this paper speaks to whether the length gain persists or whether the plate is being spent |
| direction of the claim | Sustained mTORC1 hyperactivation ARRESTS bone growth - the opposite sign to the acute bafilomycin effect. Acute and chronic lysosomal impairment are therefore not interchangeable. |

## On MTORC1 necessity

Torin1 **attenuated** and **significantly diminished** the growth effect of Baf and CQ (n=5-11 animals). It did not abolish it, and the experiment was not designed to test necessity. Three things argue against reading this as proof of mechanism:

1. Torin1 suppresses bone growth on its own (rapamycin does, and the paper cites this), so    partial attenuation is equally consistent with two independent, opposing effects.
2. **p-MTOR itself was not significantly changed** (2.0 +/- 1.2 fold, p=0.49), and    **p-RPS6KB1 was not significantly changed** in the cell model (-14 +/- 44%, p=0.78).    The only strong MTORC1 readout is p-RPS6.
3. The authors themselves flag a dose mismatch: Baf activates RPS6 ~5x more strongly than    CQ but grows bone only 24% more, which they say *'suggests that different mechanisms of    growth may be involved'*.

And their own caveat: *'we think that genetic studies are required to confirm this assumption and to extend it into a physiological setting.'* This audit therefore records MTORC1 as **contributory but not demonstrated necessary**.

## The gap that matters most

**There is no washout or recovery experiment.** The strings `washout` and `recover` do not appear anywhere in the full text. The culture ran 5-6 days with continuous exposure. Nothing in this paper addresses whether the length gain persists, whether lysosomal function recovers, or whether the plate is simply being spent faster. Every durability claim about this mechanism is currently unsupported in either direction.

## The chronic counterpart

PMID 28872463 (*J Clin Invest* 2017, 'mTORC1 hyperactivation arrests bone growth in lysosomal storage disorders by suppressing autophagy') is **publisher-restricted**: the PMC record carries only front matter and abstract, so it cannot be quantitatively extracted here and is used qualitatively only. Its direction is the opposite of the acute result - sustained mTORC1 hyperactivation *arrests* growth. Acute and chronic lysosomal impairment cannot be treated as the same intervention.
