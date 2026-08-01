# Primary screen protocol

## Assay

Normal postnatal mouse metatarsal organ culture. Metatarsals 2-4 are dissected from each hind paw, one bone per well, cultured in defined medium. The screen is **target-agnostic**: no compound is included because a pathway argument says it should work, and no compound is excluded from analysis because it lacks one.

**Primary readout: absolute longitudinal length gain, measured daily.** Not percentage change, not an endpoint measurement, not a marker.

## Controls

| control | role | concentration | basis | wells | purpose |
|---|---|---|---|---:|---|
| vehicle | VEHICLE | 0.1% DMSO, matched to the highest compound vehicle load | vehicle load is fixed across the plate, not per compound | 8 | defines the plate baseline and the smallest detectable change |
| IGF1 | PRODUCTIVE BENCHMARK | 100 ng/ml | published ex vivo metatarsal concentration (PMID 26259639) | 6 | the state-A reference: length gain without a cellular cost |
| bafilomycin A1 | TRADE-OFF BENCHMARK | 8 nM | published ex vivo metatarsal concentration (PMID 26259639) | 6 | the false-positive control - it lengthens while reducing proliferation and raising apoptosis (stage 29). Any compound whose endpoint profile matches this one has failed. |
| CNP / FGFR3-pathway stimulus | ASSAY-SENSITIVITY CONTROL | RANGE_FINDING_REQUIRED | the specific agent is chosen from the stage-49 control set; its ex vivo concentration is established on the range-finding plate before the screen runs | 6 | proves the assay can detect a growth increase at all; excluded from novelty ranking |
| antiproliferative control | CYTOTOXIC CONTROL | RANGE_FINDING_REQUIRED | concentration set to the lowest that reduces EdU index without gross explant death | 4 | defines what a proliferation-cost phenotype looks like in this assay |
| washout-only | WASHOUT CONTROL | vehicle throughout, medium changed on the washout schedule | no compound | 4 | separates the effect of the medium change and handling from the effect of withdrawal |

The bafilomycin arm is the one that makes this screen different from a length assay. It is included specifically so that a compound producing the bafilomycin endpoint profile - more length, less proliferation, more apoptosis - is recognised as a failure rather than a hit. Stage 29 showed how easy that mistake is to make from the literature alone.

## Exposure arms

| arm | schedule | why |
|---|---|---|
| continuous | compound present for the whole culture period | the default; detects sustained effects |
| short pulse | compound present for the first 48 h, then vehicle | detects transient anabolic effects that continuous exposure would mask by toxicity |
| washout + recovery | compound for the first half, then vehicle to growth cessation | the durability arm; mandatory before any compound is called a hit |

### Why the primary screen runs one arm, not three

96 compounds x 3 arms x 6 biological replicates is 1,728 treatment wells, and at six metatarsals per animal that is roughly 288 animals. That is not a pilot. The primary screen therefore runs the **continuous arm only** at 112 animals, and the pulse and washout arms are applied to Tier-1 hits in the stage-53 secondary panel - which is where the Tier-4 washout requirement sits in any case.

This is a real constraint, not a simplification for presentation. It has a cost: a compound whose only productive effect is transient will look inert in the continuous arm and will never reach the pulse arm. That failure mode is accepted explicitly, and stage 56 lists it as a reason the pilot could return a false negative.

The washout arm remains mandatory before any compound is called a hit. It has moved later in the sequence; it has not become optional.

## Design principles, and how each is implemented

| principle | implementation |
|---|---|
| randomise bones across plates | conditions and bones are shuffled with a fixed seed (20260801) before assignment, so the map is random but reproducible |
| balance litter and animal across conditions | six litters x eight animals; the six metatarsals from one animal are shuffled into different conditions and plates |
| record bone and animal identity | every well carries `animal_id`, `litter_id`, `bone_id` and `bone_position` |
| multiple bones per animal are not independent | the analysis nests bone within animal within litter; the animal is the replicate |
| plate-position controls | the outer ring of every plate is reserved for vehicle and flagged `is_edge_well`; evaporation and thermal gradient are tested as a fixed effect before any compound is read |
| blinded image analysis | wells are analysed under a scrambled identifier; the stage-51 pipeline stores the blinded key and an audit trail of every manual correction |
| technical versus biological replication | repeated daily images of one bone are technical; six bones from six different animals are the biological replicates. The model treats them as such |
| predefined exclusion rules | written below, before any data exist |

## Exclusion rules, defined in advance

An explant is excluded if, and only if:

1. it is visibly damaged at dissection (fractured, perichondrium stripped, cartilage torn);
2. its day-0 length is outside 2 SD of the litter mean, which flags a developmentally abnormal or mis-identified bone;
3. it fails to grow at all in the first 48 h in the vehicle arm, indicating a failed explant rather than a treatment effect;
4. its images fail the stage-51 quality-control flags on more than two consecutive days;
5. it is contaminated.

Exclusions are recorded per bone with the rule number and are reported in the results, including for arms where excluding helps the compound.

## Concentrations

96 of 96 pilot compounds have retrievable primary potency and get a 3x-30x range-finding bracket around it. The remaining 0 get a half-log ladder from the solubility limit. **Every compound goes through range-finding before the primary screen**, and the screen concentration is defined by a rule - the highest concentration that leaves EdU index and viability indistinguishable from vehicle - not by a number chosen now.

The vehicle load is fixed at 0.1% DMSO across the whole plate rather than varied per compound, so the top testable concentration for each compound is set by its stock solubility at 1000x. Compounds that cannot reach their bracket at that vehicle load are recorded as solubility-limited rather than tested at a higher DMSO concentration.

## What this protocol cannot deliver

Explants are avascular, unloaded, and endocrine-free. They cannot report vascular invasion, mechanical loading effects, or systemic exposure, and they cannot measure adult bone length. Nothing here supports a claim about human height, and no dosing or self-experimentation guidance appears in this or any other stage.
