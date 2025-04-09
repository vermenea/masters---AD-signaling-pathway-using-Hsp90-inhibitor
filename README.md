# Influence of 17-AAG, an Hsp90 Inhibitor, on Signaling Pathways in Atopic Dermatitis

## Introduction

The aim of this study is to analyze the impact of the Hsp90 inhibitor, 17-AAG, on cellular signaling in the context of atopic dermatitis therapy. The research focuses on understanding how 17-AAG influences key signaling proteins involved in the pathogenesis of atopic dermatitis.

## Methodology

### Cell Line and Stimulation

- **Immortalized Keratinocyte Cell Line (HaCaT):** The HaCaT cell line was used as a model for keratinocytes.
- **Stimuli:** HaCaT cells were stimulated with **IFN-γ/TNF-α** to simulate inflammatory conditions relevant to atopic dermatitis.

### Treatment Conditions

- Cells were cultured in the presence or absence of **17-AAG**, an Hsp90 inhibitor.

### Analysis

The phosphorylation status of the following signaling proteins was analyzed:
- **STAT-1** ✅
- **STAT-3** ❌todo
- **STAT-6** ✅
- **ERK**    ❌todo
- **MAPK**   ❌todo

The activation of these proteins, as well as the role of Hsp90 acetylation (acLys284/292), was assessed in activated HaCaT cells using **immunoblotting**.

Densynometry was conducted using ImageJ. The results were then imported into Excel, where they were normalized against β-actin. Subsequently, the data was transferred into Visual Studio Code for further analysis using Python.

## Results and Conclusion

Our data suggest that Hsp90 inhibition via 17-AAG has differential effects across signaling pathways:

-STAT3 and MAPK/ERK are suppressed by both 0.1 µM and 1 µM, consistent with reduced pro-inflammatory signaling and cell growth.
-STAT6, in contrast, is only upregulated at higher concentrations, possibly indicating a compensatory pathway activation or cell-specific signaling rewiring.
-These effects align with the role of Hsp90 in stabilizing a broad range of signaling proteins, including kinases and transcription factors. Its inhibition thus causes a domino-like collapse of multiple intracellular signaling cascades.
-Such dose-dependent and pathway-specific modulation makes 17-AAG a promising molecule in the context of targeted therapeutic strategies for neurodegenerative conditions, where selective suppression of inflammatory or proliferative pathways is desirable.
