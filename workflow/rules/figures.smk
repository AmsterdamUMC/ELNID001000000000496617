# Rules

rule figure1:
  input:
    sample_metadata_xlsx=config["data"]["sample_metadata_xlsx"],
  output:
    figB_pdf="output/figures/figure1/figureB.pdf",
    figC_pdf="output/figures/figure1/figureC.pdf",
  conda:
    "../envs/r-figure1.yaml"
  message:
    "--- Preparing constituent panels figure 1 ---"
  threads: 1
  resources:
    mem_mb=16000,
  log:
    "output/figures/figure1.log"
  params:
    plotting_parameters=config["R"]["plotting_parameters"],
  script:
    "../scripts/figures/figure1.R"

rule figure2:
  input:
    sample_metadata_xlsx=config["data"]["sample_metadata_xlsx"],
    dds_rds="output/rnaseq/dds/dds.Rds",
    degs_csv="output/rnaseq/analyses/BINF_RvNR_W0/degs/degs_BINF_RvNR_W0.csv",
    fgsea_csv="output/rnaseq/analyses/BINF_RvNR_W0/fgsea/fgsea_BINF_RvNR_W0.csv",
  output:
    figA_pdf="output/figures/figure2/figureA.pdf",
    figB_pdf="output/figures/figure2/figureB.pdf",
    figC_pdf="output/figures/figure2/figureC.pdf",
    figD_pdf="output/figures/figure2/figureD.pdf",
  conda:
    "../envs/r-figure2.yaml"
  message:
    "--- Preparing constituent panels figure 2 ---"
  threads: 1
  resources:
    mem_mb=16000,
  log:
    "output/figures/figure2.log"
  params:
    goi_xlsx=config["R"]["genes_xlsx"],
    pwoi_xlsx=config["R"]["genesets_xlsx"],
    plotting_parameters=config["R"]["plotting_parameters"],
  script:
    "../scripts/figures/figure2.R"

rule figure3:
  input:
    sample_metadata_xlsx=config["data"]["sample_metadata_xlsx"],
    rnaseq_dds_rds="output/rnaseq/dds/dds.Rds",
    rnaseq_vst_RvNR_W0_rds="output/rnaseq/analyses/BINF_RvNR_W0/dds/vst_BINF_RvNR_W0.Rds",
    degs_RvNR_W0_csv="output/rnaseq/analyses/BINF_RvNR_W0/degs/degs_BINF_RvNR_W0.csv",
    degs_R_W8vW0_csv="output/rnaseq/analyses/BINF_R_W8vW0/degs/degs_BINF_R_W8vW0.csv",
    degs_NR_W8vW0_csv="output/rnaseq/analyses/BINF_NR_W8vW0/degs/degs_BINF_NR_W8vW0.csv",
    fgsea_RvNR_W0_csv="output/rnaseq/analyses/BINF_RvNR_W0/fgsea/fgsea_BINF_RvNR_W0.csv",
    fgsea_R_W8vW0_csv="output/rnaseq/analyses/BINF_R_W8vW0/fgsea/fgsea_BINF_R_W8vW0.csv",
    fgsea_NR_W8vW0_csv="output/rnaseq/analyses/BINF_NR_W8vW0/fgsea/fgsea_BINF_NR_W8vW0.csv",
    olink_se_rds="output/olink/se/se.Rds",
    deps_R_W8vW0_csv="output/olink/analyses/BINF_R_W8vW0/deps/deps_BINF_R_W8vW0.csv",
    deps_NR_W8vW0_csv="output/olink/analyses/BINF_NR_W8vW0/deps/deps_BINF_NR_W8vW0.csv",
  output:
    figA_pdf="output/figures/figure3/figureA.pdf",
    figB_pdf="output/figures/figure3/figureB.pdf",
    figC_pdf="output/figures/figure3/figureC.pdf",
    figD_pdf="output/figures/figure3/figureD.pdf",
    figEG_pdf="output/figures/figure3/figureEG.pdf",
    figF_pdf="output/figures/figure3/figureF.pdf",
    figH_pdf="output/figures/figure3/figureH.pdf",
    figI_pdf="output/figures/figure3/figureI.pdf",
  conda:
    "../envs/r-figure3.yaml"
  message:
    "--- Preparing constituent panels figure 3 ---"
  threads: 1
  resources:
    mem_mb=16000,
  log:
    "output/figures/figure3.log"
  params:
    goi_xlsx=config["R"]["genes_xlsx"],
    pwoi_xlsx=config["R"]["genesets_xlsx"],
    plotting_parameters=config["R"]["plotting_parameters"],
  script:
    "../scripts/figures/figure3.R"

rule figure4:
  input:
    sample_metadata_xlsx=config["data"]["sample_metadata_xlsx"],
    vst_R_W8vW0_TOF_rds="output/rnaseq/analyses/BINF_R_W8vW0/dds/vst_BINF_R_W8vW0.Rds",
    microarray_toedter2011_R_W8vW0_IFX_se_rds="output/microarray_toedter2011/analyses/BINF_R_W8vW0/se/se_BINF_R_W8vW0.Rds",
    degs_R_W8vW0_TOF_csv="output/rnaseq/analyses/BINF_R_W8vW0/degs/degs_BINF_R_W8vW0.csv",
    fgsea_R_W8vW0_TOF_csv="output/rnaseq/analyses/BINF_R_W8vW0/fgsea/fgsea_BINF_R_W8vW0.csv",
    fgsea_NR_W8vW0_TOF_csv="output/rnaseq/analyses/BINF_NR_W8vW0/fgsea/fgsea_BINF_NR_W8vW0.csv",
    degs_R_W8vW0_TOF_melonardanaz2025_csv="output/scrnaseq_melonardanaz2025/analyses/BINF_R_W8vW0/degs/degs_BINF_R_W8vW0.csv",
    fgsea_R_W8vW0_TOF_melonardanaz2025_csv="output/scrnaseq_melonardanaz2025/analyses/BINF_R_W8vW0/fgsea/fgsea_BINF_R_W8vW0.csv",
    fgsea_NR_W8vW0_TOF_melonardanaz2025_csv="output/scrnaseq_melonardanaz2025/analyses/BINF_NR_W8vW0/fgsea/fgsea_BINF_NR_W8vW0.csv",
    degs_R_W8vW0_IFX_toedter2011_csv="output/microarray_toedter2011/analyses/BINF_R_W8vW0/degs/degs_BINF_R_W8vW0.csv",
    fgsea_R_W8vW0_IFX_toedter2011_csv="output/microarray_toedter2011/analyses/BINF_R_W8vW0/fgsea/fgsea_BINF_R_W8vW0.csv",
    fgsea_NR_W8vW0_IFX_toedter2011_csv="output/microarray_toedter2011/analyses/BINF_NR_W8vW0/fgsea/fgsea_BINF_NR_W8vW0.csv",
  output:
    figA_pdf="output/figures/figure4/figureA.pdf",
    figBD_pdf="output/figures/figure4/figureBD.pdf",
    figC_pdf="output/figures/figure4/figureC.pdf",
    figE_pdf="output/figures/figure4/figureE.pdf",
  conda:
    "../envs/r-figure4.yaml"
  message:
    "--- Preparing constituent panels figure 4 ---"
  threads: 1
  resources:
    mem_mb=16000,
  log:
    "output/figures/figure4.log"
  params:
    goi_xlsx=config["R"]["genes_xlsx"],
    pwoi_xlsx=config["R"]["genesets_xlsx"],
    plotting_parameters=config["R"]["plotting_parameters"],
  script:
    "../scripts/figures/figure4.R"

rule figure5:
  input:
    rnaseq_dds_rds="output/rnaseq/dds/dds.Rds",
    degs_BINF_R_W8vW0_csv="output/rnaseq/analyses/BINF_R_W8vW0/degs/degs_BINF_R_W8vW0.csv",
    degs_BINF_NR_W8vW0_csv="output/rnaseq/analyses/BINF_NR_W8vW0/degs/degs_BINF_NR_W8vW0.csv",
    fgsea_BINF_R_W8vW0_csv="output/rnaseq/analyses/BINF_R_W8vW0/fgsea/fgsea_BINF_R_W8vW0.csv",
    fgsea_BINF_NR_W8vW0_csv="output/rnaseq/analyses/BINF_NR_W8vW0/fgsea/fgsea_BINF_NR_W8vW0.csv",
    degs_PB_R_W8vW0_csv="output/rnaseq/analyses/PB_R_W8vW0/degs/degs_PB_R_W8vW0.csv",
    degs_PB_NR_W8vW0_csv="output/rnaseq/analyses/PB_NR_W8vW0/degs/degs_PB_NR_W8vW0.csv",
    fgsea_PB_R_W8vW0_csv="output/rnaseq/analyses/PB_R_W8vW0/fgsea/fgsea_PB_R_W8vW0.csv",
    fgsea_PB_NR_W8vW0_csv="output/rnaseq/analyses/PB_NR_W8vW0/fgsea/fgsea_PB_NR_W8vW0.csv",
  output:
    figA_pdf="output/figures/figure5/figureA.pdf",
    figB_pdf="output/figures/figure5/figureB.pdf",
    figC_pdf="output/figures/figure5/figureC.pdf",
  conda:
    "../envs/r-figure5.yaml"
  message:
    "--- Preparing constituent panels figure 5 ---"
  threads: 1
  resources:
    mem_mb=16000,
  log:
    "output/figures/figure5.log"
  params:
    goi_xlsx=config["R"]["genes_xlsx"],
    pwoi_xlsx=config["R"]["genesets_xlsx"],
    plotting_parameters=config["R"]["plotting_parameters"],
  script:
    "../scripts/figures/figure5.R"

rule supplementaryfigure1:
  input:
    deps_R_W8vW0_csv="output/olink/analyses/BINF_R_W8vW0/deps/deps_BINF_R_W8vW0.csv",
    deps_NR_W8vW0_csv="output/olink/analyses/BINF_NR_W8vW0/deps/deps_BINF_NR_W8vW0.csv",
  output:
    figA_pdf="output/figures/supplementary_figure1/figureA.pdf",
    figB_pdf="output/figures/supplementary_figure1/figureB.pdf",
    figC_pdf="output/figures/supplementary_figure1/figureC.pdf",
  conda:
    "../envs/r-supplementary_figure1.yaml"
  message:
    "--- Preparing constituent panels supplementary figure 1 ---"
  threads: 1
  resources:
    mem_mb=16000,
  log:
    "output/figures/supplementaryfigure1.log"
  params:
    plotting_parameters=config["R"]["plotting_parameters"],
  script:
    "../scripts/figures/supplementary_figure1.R"

rule supplementaryfigure2:
  input:
    degs_NR_W8vW0_TOF_csv="output/rnaseq/analyses/BINF_NR_W8vW0/degs/degs_BINF_NR_W8vW0.csv",
    fgsea_R_W8vW0_TOF_csv="output/rnaseq/analyses/BINF_R_W8vW0/fgsea/fgsea_BINF_R_W8vW0.csv",
    fgsea_NR_W8vW0_TOF_csv="output/rnaseq/analyses/BINF_NR_W8vW0/fgsea/fgsea_BINF_NR_W8vW0.csv",
    degs_NR_W8vW0_TOF_melonardanaz2025_csv="output/scrnaseq_melonardanaz2025/analyses/BINF_NR_W8vW0/degs/degs_BINF_NR_W8vW0.csv",
    fgsea_R_W8vW0_TOF_melonardanaz2025_csv="output/scrnaseq_melonardanaz2025/analyses/BINF_R_W8vW0/fgsea/fgsea_BINF_R_W8vW0.csv",
    fgsea_NR_W8vW0_TOF_melonardanaz2025_csv="output/scrnaseq_melonardanaz2025/analyses/BINF_NR_W8vW0/fgsea/fgsea_BINF_NR_W8vW0.csv",
  output:
    figA_pdf="output/figures/supplementary_figure2/figureA.pdf",
    figB_pdf="output/figures/supplementary_figure2/figureB.pdf",
    figC_pdf="output/figures/supplementary_figure2/figureC.pdf",
  conda:
    "../envs/r-supplementary_figure2.yaml"
  message:
    "--- Preparing constituent panels supplementary figure 2 ---"
  threads: 1
  resources:
    mem_mb=16000,
  log:
    "output/figures/supplementaryfigure2.log"
  params:
    plotting_parameters=config["R"]["plotting_parameters"],
  script:
    "../scripts/figures/supplementary_figure2.R"
