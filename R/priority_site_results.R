################################ Priority sites ################################
#species + weights----
#Get cell values for Kolmogorov-smirnoff tests
PM_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, PM_spec_bin)
PM_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, PM_KBA_bin)
PM_RAN_vals <- get_msk_vals(RAN_var_ras, PM_RAN_bin)

VG_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, VG_spec_bin)
VG_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, VG_KBA_bin)
VG_RAN_vals <- get_msk_vals(RAN_var_ras, VG_RAN_bin)

DG_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, DG_spec_bin)
DG_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, DG_KBA_bin)
DG_RAN_vals <- get_msk_vals(RAN_var_ras, DG_RAN_bin)

TB_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, TB_spec_bin)
TB_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, TB_KBA_bin)
TB_RAN_vals <- get_msk_vals(RAN_var_ras, TB_RAN_bin)

PL_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, PL_spec_bin)
PL_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, PL_KBA_bin)
PL_RAN_vals <- get_msk_vals(RAN_var_ras, PL_RAN_bin)

AM_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, AM_spec_bin)
AM_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, AM_KBA_bin)
AM_RAN_vals <- get_msk_vals(RAN_var_ras, AM_RAN_bin)

MF_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, MF_spec_bin)
MF_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, MF_KBA_bin)
MF_RAN_vals <- get_msk_vals(RAN_var_ras, MF_RAN_bin)

MD_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, MD_spec_bin)
MD_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, MD_KBA_bin)
MD_RAN_vals <- get_msk_vals(RAN_var_ras, MD_RAN_bin)

LH_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, LH_spec_bin)
LH_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, LH_KBA_bin)
LH_RAN_vals <- get_msk_vals(RAN_var_ras, LH_RAN_bin)

VV_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, VV_spec_bin)
VV_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, VV_KBA_bin)
VV_RAN_vals <- get_msk_vals(RAN_var_ras, VV_RAN_bin)

BT_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, BT_spec_bin)
BT_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, BT_KBA_bin)
BT_RAN_vals <- get_msk_vals(RAN_var_ras, BT_RAN_bin)

HA_spec_vals <- get_msk_vals(CAZ_wgt_var_ras, HA_spec_bin)
HA_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, HA_KBA_bin)
HA_RAN_vals <- get_msk_vals(RAN_var_ras, HA_RAN_bin)

#Create violin plots
PM_vals <- c(PM_spec_vals, PM_KBA_vals)
VG_vals <- c(VG_spec_vals, VG_KBA_vals)
DG_vals <- c(DG_spec_vals, DG_KBA_vals)
TB_vals <- c(TB_spec_vals, TB_KBA_vals)
PL_vals <- c(PL_spec_vals, PL_KBA_vals)
AM_vals <- c(AM_spec_vals, AM_KBA_vals)
MF_vals <- c(MF_spec_vals, MF_KBA_vals)
MD_vals <- c(MD_spec_vals, MD_KBA_vals)
LH_vals <- c(LH_spec_vals, LH_KBA_vals)
VV_vals <- c(VV_spec_vals, VV_KBA_vals)
BT_vals <- c(BT_spec_vals, BT_KBA_vals)
HA_vals <- c(HA_spec_vals, HA_KBA_vals)
vals <- c(PM_vals, VG_vals,DG_vals,TB_vals,PL_vals,AM_vals,MF_vals,MD_vals,LH_vals,VV_vals,BT_vals, HA_vals)
nms <- c(rep("Pheidole megacephala", length(PM_vals)),
         rep("Vespula germanica", length(VG_vals)),
         rep("Digitonthophagus gazella", length(DG_vals)),
         rep("Tetramorium bicarinatum", length(TB_vals)),
         rep("Paratrechina longicornis", length(PL_vals)),
         rep("Apis mellifera", length(AM_vals)),
         rep("Monomorium floricola", length(MF_vals)),
         rep("Monomorium destructor", length(MD_vals)),
         rep("Linepithema humile", length(LH_vals)),
         rep("Vespula vulgaris", length(VV_vals)),
         rep("Bombus terrestris", length(BT_vals)),
         rep("Heteronychus arator", length(HA_vals)))
type <- c(rep("species", length(PM_spec_vals)), 
          rep("KBA", length(PM_KBA_vals)),
          rep("species", length(VG_spec_vals)), 
          rep("KBA", length(VG_KBA_vals)),
          rep("species", length(DG_spec_vals)), 
          rep("KBA", length(DG_KBA_vals)),
          rep("species", length(TB_spec_vals)), 
          rep("KBA", length(TB_KBA_vals)),
          rep("species", length(PL_spec_vals)), 
          rep("KBA", length(PL_KBA_vals)),
          rep("species", length(AM_spec_vals)), 
          rep("KBA", length(AM_KBA_vals)),
          rep("species", length(MF_spec_vals)), 
          rep("KBA", length(MF_KBA_vals)),
          rep("species", length(MD_spec_vals)), 
          rep("KBA", length(MD_KBA_vals)),
          rep("species", length(LH_spec_vals)), 
          rep("KBA", length(LH_KBA_vals)),
          rep("species", length(VV_spec_vals)), 
          rep("KBA", length(VV_KBA_vals)),
          rep("species", length(BT_spec_vals)), 
          rep("KBA", length(BT_KBA_vals)),
          rep("species", length(HA_spec_vals)), 
          rep("KBA", length(HA_KBA_vals)))
meds <- c(rep(median(PM_spec_vals), length(PM_spec_vals)), rep(median(PM_KBA_vals), length(PM_KBA_vals)),
          rep(median(VG_spec_vals), length(VG_spec_vals)), rep(median(VG_KBA_vals), length(VG_KBA_vals)),
          rep(median(DG_spec_vals), length(DG_spec_vals)), rep(median(DG_KBA_vals), length(DG_KBA_vals)),
          rep(median(TB_spec_vals), length(TB_spec_vals)), rep(median(TB_KBA_vals), length(TB_KBA_vals)),
          rep(median(PL_spec_vals), length(PL_spec_vals)), rep(median(PL_KBA_vals), length(PL_KBA_vals)),
          rep(median(AM_spec_vals), length(AM_spec_vals)), rep(median(AM_KBA_vals), length(AM_KBA_vals)),
          rep(median(MF_spec_vals), length(MF_spec_vals)), rep(median(MF_KBA_vals), length(MF_KBA_vals)),
          rep(median(MD_spec_vals), length(MD_spec_vals)), rep(median(MD_KBA_vals), length(MD_KBA_vals)),
          rep(median(LH_spec_vals), length(LH_spec_vals)), rep(median(LH_KBA_vals), length(LH_KBA_vals)),
          rep(median(VV_spec_vals), length(VV_spec_vals)), rep(median(VV_KBA_vals), length(VV_KBA_vals)),
          rep(median(BT_spec_vals), length(BT_spec_vals)), rep(median(BT_KBA_vals), length(BT_KBA_vals)),
          rep(median(HA_spec_vals), length(HA_spec_vals)), rep(median(HA_KBA_vals), length(HA_KBA_vals)))
df <- data.frame(nms, type, vals, meds)

val_plot <- ggplot(df, aes(x = vals, y = nms, fill = type)) +
  geom_violin(draw_quantiles = 0.5, 
              adjust = 0.2, #adjust changes the smoothness; lower is more faithful to the data
              scale = "width") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.text.y = element_text(face = "italic", size = 12),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  scale_fill_manual(values = c("#66B2FF", "#FFB266"),
                    labels = c("KBA mask", "species + weight"),
                    name = "Scenarios") +
  scale_y_discrete(limits = rev) +
  xlab("Site sensitivty")


#Kolmogorov-smirnoff tests - KBA vs no KBA
#a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution
#presence of ties a result of rounding
PM_KBA <- ks.test(PM_KBA_vals, PM_spec_vals, alternative = "two.sided")
VG_KBA <- ks.test(VG_KBA_vals, VG_spec_vals, alternative = "two.sided")
DG_KBA <- ks.test(DG_KBA_vals, DG_spec_vals, alternative = "two.sided")
TB_KBA <- ks.test(TB_KBA_vals, TB_spec_vals, alternative = "two.sided")
PL_KBA <- ks.test(PL_KBA_vals, PL_spec_vals, alternative = "two.sided")
AM_KBA <- ks.test(AM_KBA_vals, AM_spec_vals, alternative = "two.sided")
MF_KBA <- ks.test(MF_KBA_vals, MF_spec_vals, alternative = "two.sided")
MD_KBA <- ks.test(MD_KBA_vals, MD_spec_vals, alternative = "two.sided")
LH_KBA <- ks.test(LH_KBA_vals, LH_spec_vals, alternative = "two.sided")
VV_KBA <- ks.test(VV_KBA_vals, VV_spec_vals, alternative = "two.sided")
BT_KBA <- ks.test(BT_KBA_vals, BT_spec_vals, alternative = "two.sided")
HA_KBA <- ks.test(HA_KBA_vals, HA_spec_vals, alternative = "two.sided")

kbas <- tribble(
  ~species, ~statistic,
  "Pheidole megacephala", paste0(round(PM_KBA$statistic,4),"***"), 
  "Vespula germanica", paste0(round(VG_KBA$statistic,4),"***"), 
  "Digitonthophagus gazella", paste0(round(DG_KBA$statistic,4),"***"), 
  "Paratrechina longicornis", paste0(round(PL_KBA$statistic,4),"***"), 
  "Tetramorium bicarinatum", paste0(round(TB_KBA$statistic,4),"***"), 
  "Apis mellifera", paste0(round(AM_KBA$statistic,4),"***"), 
  "Monomorium floricola", paste0(round(MF_KBA$statistic,4),"***"), 
  "Monomorium destructor", paste0(round(MD_KBA$statistic,4),"***"), 
  "Linepithema humile", paste0(round(LH_KBA$statistic,4),"***"),
  "Vespula vulgaris", paste0(round(VV_KBA$statistic,4),"***"),
  "Bombus terrestris", paste0(round(BT_KBA$statistic,4),"***"), 
  "Heteronychus arator", paste0(round(HA_KBA$statistic,4),"***"), 
)

#Difference between IAS - species + weight
PM_VG <- ks.test(PM_spec_vals, VG_spec_vals, alternative = "two.sided")
PM_DG <- ks.test(PM_spec_vals, DG_spec_vals, alternative = "two.sided")
PM_TB <- ks.test(PM_spec_vals, TB_spec_vals, alternative = "two.sided")
PM_PL <- ks.test(PM_spec_vals, PL_spec_vals, alternative = "two.sided")
PM_AM <- ks.test(PM_spec_vals, AM_spec_vals, alternative = "two.sided")
PM_MF <- ks.test(PM_spec_vals, MF_spec_vals, alternative = "two.sided")
PM_MD <- ks.test(PM_spec_vals, MD_spec_vals, alternative = "two.sided")
PM_LH <- ks.test(PM_spec_vals, LH_spec_vals, alternative = "two.sided")
PM_VV <- ks.test(PM_spec_vals, VV_spec_vals, alternative = "two.sided")
PM_BT <- ks.test(PM_spec_vals, BT_spec_vals, alternative = "two.sided")
PM_HA <- ks.test(PM_spec_vals, HA_spec_vals, alternative = "two.sided")

VG_DG <- ks.test(VG_spec_vals, DG_spec_vals, alternative = "two.sided")
VG_TB <- ks.test(VG_spec_vals, TB_spec_vals, alternative = "two.sided")
VG_PL <- ks.test(VG_spec_vals, PL_spec_vals, alternative = "two.sided")
VG_AM <- ks.test(VG_spec_vals, AM_spec_vals, alternative = "two.sided")
VG_MF <- ks.test(VG_spec_vals, MF_spec_vals, alternative = "two.sided")
VG_MD <- ks.test(VG_spec_vals, MD_spec_vals, alternative = "two.sided")
VG_LH <- ks.test(VG_spec_vals, LH_spec_vals, alternative = "two.sided")
VG_VV <- ks.test(VG_spec_vals, VV_spec_vals, alternative = "two.sided")
VG_BT <- ks.test(VG_spec_vals, BT_spec_vals, alternative = "two.sided")
VG_HA <- ks.test(VG_spec_vals, HA_spec_vals, alternative = "two.sided")

DG_TB <- ks.test(DG_spec_vals, TB_spec_vals, alternative = "two.sided")
DG_PL <- ks.test(DG_spec_vals, PL_spec_vals, alternative = "two.sided")
DG_AM <- ks.test(DG_spec_vals, AM_spec_vals, alternative = "two.sided")
DG_MF <- ks.test(DG_spec_vals, MF_spec_vals, alternative = "two.sided")
DG_MD <- ks.test(DG_spec_vals, MD_spec_vals, alternative = "two.sided")
DG_LH <- ks.test(DG_spec_vals, LH_spec_vals, alternative = "two.sided")
DG_VV <- ks.test(DG_spec_vals, VV_spec_vals, alternative = "two.sided")
DG_BT <- ks.test(DG_spec_vals, BT_spec_vals, alternative = "two.sided")
DG_HA <- ks.test(DG_spec_vals, HA_spec_vals, alternative = "two.sided")

TB_PL <- ks.test(TB_spec_vals, PL_spec_vals, alternative = "two.sided")
TB_AM <- ks.test(TB_spec_vals, AM_spec_vals, alternative = "two.sided")
TB_MF <- ks.test(TB_spec_vals, MF_spec_vals, alternative = "two.sided")
TB_MD <- ks.test(TB_spec_vals, MD_spec_vals, alternative = "two.sided")
TB_LH <- ks.test(TB_spec_vals, LH_spec_vals, alternative = "two.sided")
TB_VV <- ks.test(TB_spec_vals, VV_spec_vals, alternative = "two.sided")
TB_BT <- ks.test(TB_spec_vals, BT_spec_vals, alternative = "two.sided")
TB_HA <- ks.test(TB_spec_vals, HA_spec_vals, alternative = "two.sided")

PL_AM <- ks.test(PL_spec_vals, AM_spec_vals, alternative = "two.sided")
PL_MF <- ks.test(PL_spec_vals, MF_spec_vals, alternative = "two.sided")
PL_MD <- ks.test(PL_spec_vals, MD_spec_vals, alternative = "two.sided")
PL_LH <- ks.test(PL_spec_vals, LH_spec_vals, alternative = "two.sided")
PL_VV <- ks.test(PL_spec_vals, VV_spec_vals, alternative = "two.sided")
PL_BT <- ks.test(PL_spec_vals, BT_spec_vals, alternative = "two.sided")
PL_HA <- ks.test(PL_spec_vals, HA_spec_vals, alternative = "two.sided")

AM_MF <- ks.test(AM_spec_vals, MF_spec_vals, alternative = "two.sided")
AM_MD <- ks.test(AM_spec_vals, MD_spec_vals, alternative = "two.sided")
AM_LH <- ks.test(AM_spec_vals, LH_spec_vals, alternative = "two.sided")
AM_VV <- ks.test(AM_spec_vals, VV_spec_vals, alternative = "two.sided")
AM_BT <- ks.test(AM_spec_vals, BT_spec_vals, alternative = "two.sided")
AM_HA <- ks.test(AM_spec_vals, HA_spec_vals, alternative = "two.sided")

MF_MD <- ks.test(MF_spec_vals, MD_spec_vals, alternative = "two.sided")
MF_LH <- ks.test(MF_spec_vals, LH_spec_vals, alternative = "two.sided")
MF_VV <- ks.test(MF_spec_vals, VV_spec_vals, alternative = "two.sided")
MF_BT <- ks.test(MF_spec_vals, BT_spec_vals, alternative = "two.sided")
MF_HA <- ks.test(MF_spec_vals, HA_spec_vals, alternative = "two.sided")

MD_LH <- ks.test(MD_spec_vals, LH_spec_vals, alternative = "two.sided")
MD_VV <- ks.test(MD_spec_vals, VV_spec_vals, alternative = "two.sided")
MD_BT <- ks.test(MD_spec_vals, BT_spec_vals, alternative = "two.sided")
MD_HA <- ks.test(MD_spec_vals, HA_spec_vals, alternative = "two.sided")

LH_VV <- ks.test(LH_spec_vals, VV_spec_vals, alternative = "two.sided")
LH_BT <- ks.test(LH_spec_vals, BT_spec_vals, alternative = "two.sided")
LH_HA <- ks.test(LH_spec_vals, HA_spec_vals, alternative = "two.sided")

VV_BT <- ks.test(VV_spec_vals, BT_spec_vals, alternative = "two.sided")
VV_HA <- ks.test(VV_spec_vals, HA_spec_vals, alternative = "two.sided")

BT_HA <- ks.test(BT_spec_vals, HA_spec_vals, alternative = "two.sided")

v <- c(0,PM_VG$statistic,PM_DG$statistic,PM_TB$statistic,PM_PL$statistic,PM_AM$statistic,PM_MF$statistic,PM_MD$statistic,PM_LH$statistic,PM_VV$statistic,PM_BT$statistic,PM_HA$statistic,
       0,0,VG_DG$statistic,VG_TB$statistic,VG_PL$statistic,VG_AM$statistic,VG_MF$statistic,VG_MD$statistic,VG_LH$statistic,VG_VV$statistic,VG_BT$statistic,VG_HA$statistic,
       0,0,0,DG_TB$statistic,DG_PL$statistic,DG_AM$statistic,DG_MF$statistic,DG_MD$statistic,DG_LH$statistic,DG_VV$statistic,DG_BT$statistic,DG_HA$statistic,
       0,0,0,0,TB_PL$statistic,TB_AM$statistic,TB_MF$statistic,TB_MD$statistic,TB_LH$statistic,TB_VV$statistic,TB_BT$statistic,TB_HA$statistic,
       0,0,0,0,0,PL_AM$statistic,PL_MF$statistic,PL_MD$statistic,PL_LH$statistic,PL_VV$statistic,PL_BT$statistic,PL_HA$statistic,
       0,0,0,0,0,0,AM_MF$statistic,AM_MD$statistic,AM_LH$statistic,AM_VV$statistic,AM_BT$statistic,AM_HA$statistic,
       0,0,0,0,0,0,0,MF_MD$statistic,MF_LH$statistic,MF_VV$statistic,MF_BT$statistic,MF_HA$statistic,
       0,0,0,0,0,0,0,0,MD_LH$statistic,MD_VV$statistic,MD_BT$statistic,MD_HA$statistic,
       0,0,0,0,0,0,0,0,0,LH_VV$statistic,LH_BT$statistic,LH_HA$statistic,
       0,0,0,0,0,0,0,0,0,0,VV_BT$statistic,VV_HA$statistic,
       0,0,0,0,0,0,0,0,0,0,0,BT_HA$statistic,
       0,0,0,0,0,0,0,0,0,0,0,0)

tm <- matrix(v, nrow = 12, ncol = 12)
rownames(tm) <- c(":italic(Pheidole~~megacephala)", 
                  ":italic(Vespula~~germanica)", 
                  ":italic(Digitonthophagus~~gazella)", 
                  ":italic(Tetramorium~~bicarinatum)", 
                  ":italic(Paratrechina~~longicornis)",
                  ":italic(Apis~~mellifera)",
                  ":italic(Monomorium~~floricola)",
                  ":italic(Monomorium~~destructor)",
                  ":italic(Linepithema~~humile)",
                  ":italic(Vespula~~vulgaris)",
                  ":italic(Bombus~~terrestris)",
                  ":italic(Heteronychus~~arator)")
colnames(tm) <- c("Pm", "Vg", "Dg", "Tb", "Pl", "Am", "Mf", "Md", "Lh", "Vv", "Bt", "HA")
cp <- corrplot::corrplot(tm, type = "lower", method = "color", 
                   cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
                   tl.srt = 0, tl.col = "black", tl.cex = 0.8, addCoef.col = "black",
                   mar = c(0,0,0,0))
ggsave(filename = "spec_only_ks.pdf", plot = cp, device = "pdf", path = file.path("Outcome", "pdf"))

#Proportion difference (KBA vs no KBA vs Random) in number of top sensitive sites
#Top two (i.e. >= 0.98 sensitivity)
props <- c(1, 0.98, 0.95, 0.90, 0.75, 0.50, 0.25, 0.00)
diffs <- c()

PM_spec_prop <- multi_props(PM_spec_vals, props)
PM_KBA_prop <- multi_props(PM_KBA_vals, props)
PM_RAN_prop <- multi_props(PM_RAN_vals, props)

VG_spec_prop <- multi_props(VG_spec_vals, props)
VG_KBA_prop <- multi_props(VG_KBA_vals, props)
VG_RAN_prop <- multi_props(VG_RAN_vals, props)

DG_spec_prop <- multi_props(DG_spec_vals, props)
DG_KBA_prop <- multi_props(DG_KBA_vals, props)
DG_RAN_prop <- multi_props(DG_RAN_vals, props)

TB_spec_prop <- multi_props(TB_spec_vals, props)
TB_KBA_prop <- multi_props(TB_KBA_vals, props)
TB_RAN_prop <- multi_props(TB_RAN_vals, props)

PL_spec_prop <- multi_props(PL_spec_vals, props)
PL_KBA_prop <- multi_props(PL_KBA_vals, props)
PL_RAN_prop <- multi_props(PL_RAN_vals, props)

AM_spec_prop <- multi_props(AM_spec_vals, props)
AM_KBA_prop <- multi_props(AM_KBA_vals, props)
AM_RAN_prop <- multi_props(AM_RAN_vals, props)

MF_spec_prop <- multi_props(MF_spec_vals, props)
MF_KBA_prop <- multi_props(MF_KBA_vals, props)
MF_RAN_prop <- multi_props(MF_RAN_vals, props)

MD_spec_prop <- multi_props(MD_spec_vals, props)
MD_KBA_prop <- multi_props(MD_KBA_vals, props)
MD_RAN_prop <- multi_props(MD_RAN_vals, props)

LH_spec_prop <- multi_props(LH_spec_vals, props)
LH_KBA_prop <- multi_props(LH_KBA_vals, props)
LH_RAN_prop <- multi_props(LH_RAN_vals, props)

VV_spec_prop <- multi_props(VV_spec_vals, props)
VV_KBA_prop <- multi_props(VV_KBA_vals, props)
VV_RAN_prop <- multi_props(VV_RAN_vals, props)

BT_spec_prop <- multi_props(BT_spec_vals, props)
BT_KBA_prop <- multi_props(BT_KBA_vals, props)
BT_RAN_prop <- multi_props(BT_RAN_vals, props)

HA_spec_prop <- multi_props(HA_spec_vals, props)
HA_KBA_prop <- multi_props(HA_KBA_vals, props)
HA_RAN_prop <- multi_props(HA_RAN_vals, props)

Total <- as.data.frame(rbind(PM_spec_prop,VG_spec_prop,DG_spec_prop,TB_spec_prop,PL_spec_prop,AM_spec_prop,MF_spec_prop,MD_spec_prop,LH_spec_prop,VV_spec_prop,BT_spec_prop,HA_spec_prop,
                             PM_KBA_prop, VG_KBA_prop, DG_KBA_prop, TB_KBA_prop, PL_KBA_prop,AM_KBA_prop,MF_KBA_prop,MD_KBA_prop,LH_KBA_prop,VV_KBA_prop,BT_KBA_prop,HA_KBA_prop,
                             PM_RAN_prop, VG_RAN_prop, DG_RAN_prop, TB_RAN_prop, PL_RAN_prop,AM_RAN_prop,MF_RAN_prop,MD_RAN_prop,LH_RAN_prop,VV_RAN_prop,BT_RAN_prop,HA_RAN_prop))
Type <- c(rep("Unmasked", 12), rep("Masked", 12), rep("Random", 12))
Species <- rep(c("P. megacephala", "V. germanica", "D. gazella",
                 "T. bicarinatum", "P. longicornis", "A. mellifera",
                 "M. floricola", "M. destructor", "L. humile",
                 "V. vulgaris", "B. terrestris", "H. arator"), 3)
Total <- cbind(Species, Type,Total)
Total <- Total %>%
  as_tibble() %>%
  mutate(Species = factor(Species)) %>%
  mutate(Type = factor(Type)) %>%
  pivot_longer(cols = !Species & !Type, 
               names_to = "SiteSensitivity", 
               values_to = "DistributionCoverage") %>%
  mutate(SiteSensitivity = as.double(SiteSensitivity))
Total_min <- Total %>%
  dplyr::filter(SiteSensitivity == 1 | SiteSensitivity == 0.98)

df <- data.frame(xmin = c(1,0.98,0.95,0.9,0.75,0.5,0.25),
                 xmax = c(0.98,0.95,0.9,0.75,0.5,0.25,0),
                 ymin = rep(-Inf,7),
                 ymax = rep(Inf,7))

p <- ggplot(Total, 
            aes(x = SiteSensitivity, 
                y = DistributionCoverage, 
                colour = Species,
                group = interaction(Species,Type))) +
  geom_line(stat = "identity", aes(linetype = Type), size = 0.8) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) + #prob change palette
  ylab("Distribution Coverage (%)") +
  xlab("Site Sensitivity") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0))
p2 <-  p + 
  geom_rect(data = df, aes(xmin = rev(xmin), xmax = rev(xmax), ymin=ymin, ymax=ymax),
            fill = leg$colors, alpha = 0.2, inherit.aes = F) +
  font("legend.text", face = "italic")

#For inset
ins <- ggplot(Total_min, 
              aes(x = SiteSensitivity, 
                  y = DistributionCoverage, 
                  colour = Species,
                  group = interaction(Species,Type))) +
  geom_line(stat = "identity", aes(linetype = Type), size = 0.8) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) +
  ylab("Distribution Coverage (%)") +
  xlab("Site Sensitivity") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = alpha(leg$colors[7],0.2)),
        axis.text = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank()) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0.98,1)) +
  scale_y_continuous(expand = c(0,0))

#An accompanying figure could be showing numerical difference between non-kba and kba
#Take differences, e.g. Here negative values means higher distribution coverage for KBAs
pm_diff <- PM_spec_prop - PM_KBA_prop
vg_diff <- VG_spec_prop - VG_KBA_prop
dg_diff <- DG_spec_prop - DG_KBA_prop
tb_diff <- TB_spec_prop - TB_KBA_prop
pl_diff <- PL_spec_prop - PL_KBA_prop
am_diff <- AM_spec_prop - AM_KBA_prop
mf_diff <- MF_spec_prop - MF_KBA_prop
md_diff <- MD_spec_prop - MD_KBA_prop
lh_diff <- LH_spec_prop - LH_KBA_prop
vv_diff <- VV_spec_prop - VV_KBA_prop
bt_diff <- BT_spec_prop - BT_KBA_prop
ha_diff <- HA_spec_prop - HA_KBA_prop

#Then just follow through like before but using these values
Total_diff <- data.frame(rbind(pm_diff, vg_diff, dg_diff, tb_diff, pl_diff, 
                               am_diff, mf_diff, md_diff, lh_diff, vv_diff,
                               bt_diff, ha_diff))
Species <- c("Pheidole megacephala", "Vespula germanica", "Digitonthophagus gazella",
             "Tetramorium bicarinatum", "Paratrechina longicornis",
             "Apis mellifera", "Monomorium floricola", "Monomorium destructor",
             "Linepithema humile", "Vespula vulgaris", "Bombus terrestris", 
             "Heteronychus arator")
Total_diff <- data.frame(Species, Total_diff)
names(Total_diff) <- c("Species", "1.00","0.98","0.95","0.90","0.75","0.50","0.25", "0.00")
Total_diff <- Total_diff %>%
  as_tibble() %>%
  mutate(Species = factor(Species)) %>%
  pivot_longer(cols = !Species, 
               names_to = "SiteSensitivity", 
               values_to = "ScenarioDifference") %>%
  mutate(SiteSensitivity = as.double(SiteSensitivity))

Total_diff <- Total_diff %>% 
  group_by(SiteSensitivity) %>% 
  mutate(mu = mean(ScenarioDifference)) %>% 
  mutate(sd = sqrt(var(ScenarioDifference)))

ggplot(Total_diff, 
       aes(x = SiteSensitivity, 
           y = ScenarioDifference, 
           colour = Species)) +
  geom_line(stat = "identity", size = 0.8) +
  scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  geom_hline(yintercept = 0) +
  geom_rect(data = df, aes(xmin = rev(xmin), xmax = rev(xmax), ymin=ymin, ymax=ymax),
             fill = leg$colors, alpha = 0.2, inherit.aes = F) +
  font("legend.text", face = "italic") +
  ylab("Scenario difference in distribution coverage") +
  xlab("Site sensitivity")

#Mean difference among species
ggplot(Total_diff, 
       aes(x = SiteSensitivity, 
           y = mu)) +
  geom_smooth(col = "black") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  geom_hline(yintercept = 0) +
  geom_rect(data = df, aes(xmin = rev(xmin), xmax = rev(xmax), ymin=ymin, ymax=ymax),
            fill = leg$colors, alpha = 0.2, inherit.aes = F) +
  font("legend.text", face = "italic") +
  ylab("Scenario difference in distribution coverage") +
  xlab("Site sensitivity")



#species + area + weights----
#Get cell values for Kolmogorov-smirnoff tests
PM_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, PM_area_bin)
PM_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, PM_area_KBA_bin)
PM_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, PM_area_RAN_bin)

VG_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, VG_area_bin)
VG_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, VG_area_KBA_bin)
VG_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, VG_area_RAN_bin)

DG_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, DG_area_bin)
DG_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, DG_area_KBA_bin)
DG_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, DG_area_RAN_bin)

TB_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, TB_area_bin)
TB_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, TB_area_KBA_bin)
TB_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, TB_area_RAN_bin)

PL_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, PL_area_bin)
PL_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, PL_area_KBA_bin)
PL_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, PL_area_RAN_bin)

AM_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, AM_area_bin)
AM_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, AM_area_KBA_bin)
AM_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, AM_area_RAN_bin)

MF_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MF_area_bin)
MF_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MF_area_KBA_bin)
MF_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, MF_area_RAN_bin)

MD_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MD_area_bin)
MD_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MD_area_KBA_bin)
MD_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, MD_area_RAN_bin)

LH_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, LH_area_bin)
LH_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, LH_area_KBA_bin)
LH_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, LH_area_RAN_bin)

VV_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, VV_area_bin)
VV_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, VV_area_KBA_bin)
VV_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, VV_area_RAN_bin)

BT_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, BT_area_bin)
BT_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, BT_area_KBA_bin)
BT_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, BT_area_RAN_bin)

HA_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, HA_area_bin)
HA_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, HA_area_KBA_bin)
HA_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, HA_area_RAN_bin)

#Create violin plots
PM_vals <- c(PM_area_vals, PM_area_KBA_vals)
VG_vals <- c(VG_area_vals, VG_area_KBA_vals)
DG_vals <- c(DG_area_vals, DG_area_KBA_vals)
TB_vals <- c(TB_area_vals, TB_area_KBA_vals)
PL_vals <- c(PL_area_vals, PL_area_KBA_vals)
AM_vals <- c(AM_area_vals, AM_area_KBA_vals)
MF_vals <- c(MF_area_vals, MF_area_KBA_vals)
MD_vals <- c(MD_area_vals, MD_area_KBA_vals)
LH_vals <- c(LH_area_vals, LH_area_KBA_vals)
VV_vals <- c(VV_area_vals, VV_area_KBA_vals)
BT_vals <- c(BT_area_vals, BT_area_KBA_vals)
HA_vals <- c(HA_area_vals, HA_area_KBA_vals)
vals <- c(PM_vals, VG_vals,DG_vals,TB_vals,PL_vals,
          AM_vals,MF_vals,MD_vals,LH_vals,VV_vals,
          BT_vals, HA_vals)
nms <- c(rep("Pheidole megacephala", length(PM_vals)),
         rep("Vespula germanica", length(VG_vals)),
         rep("Digitonthophagus gazella", length(DG_vals)),
         rep("Tetramorium bicarinatum", length(TB_vals)),
         rep("Paratrechina longicornis", length(PL_vals)),
         rep("Apis mellifera", length(AM_vals)),
         rep("Monomorium floricola", length(MF_vals)),
         rep("Monomorium destructor", length(MD_vals)),
         rep("Linepithema humile", length(LH_vals)),
         rep("Vespula vulgaris", length(VV_vals)),
         rep("Bombus terrestris", length(BT_vals)),
         rep("Heteronychus arator", length(HA_vals)))
type <- c(rep("species", length(PM_area_vals)), 
          rep("KBA", length(PM_area_KBA_vals)),
          rep("species", length(VG_area_vals)), 
          rep("KBA", length(VG_area_KBA_vals)),
          rep("species", length(DG_area_vals)), 
          rep("KBA", length(DG_area_KBA_vals)),
          rep("species", length(TB_area_vals)), 
          rep("KBA", length(TB_area_KBA_vals)),
          rep("species", length(PL_area_vals)), 
          rep("KBA", length(PL_area_KBA_vals)),
          rep("species", length(AM_area_vals)), 
          rep("KBA", length(AM_area_KBA_vals)),
          rep("species", length(MF_area_vals)), 
          rep("KBA", length(MF_area_KBA_vals)),
          rep("species", length(MD_area_vals)), 
          rep("KBA", length(MD_area_KBA_vals)),
          rep("species", length(LH_area_vals)), 
          rep("KBA", length(LH_area_KBA_vals)),
          rep("species", length(VV_area_vals)), 
          rep("KBA", length(VV_area_KBA_vals)),
          rep("species", length(BT_area_vals)), 
          rep("KBA", length(BT_area_KBA_vals)),
          rep("species", length(HA_area_vals)), 
          rep("KBA", length(HA_area_KBA_vals)))
meds <- c(rep(median(PM_area_vals), length(PM_area_vals)), rep(median(PM_area_KBA_vals), length(PM_area_KBA_vals)),
          rep(median(VG_area_vals), length(VG_area_vals)), rep(median(VG_area_KBA_vals), length(VG_area_KBA_vals)),
          rep(median(DG_area_vals), length(DG_area_vals)), rep(median(DG_area_KBA_vals), length(DG_area_KBA_vals)),
          rep(median(TB_area_vals), length(TB_area_vals)), rep(median(TB_area_KBA_vals), length(TB_area_KBA_vals)),
          rep(median(PL_area_vals), length(PL_area_vals)), rep(median(PL_area_KBA_vals), length(PL_area_KBA_vals)),
          rep(median(AM_area_vals), length(AM_area_vals)), rep(median(AM_area_KBA_vals), length(AM_area_KBA_vals)),
          rep(median(MF_area_vals), length(MF_area_vals)), rep(median(MF_area_KBA_vals), length(MF_area_KBA_vals)),
          rep(median(MD_area_vals), length(MD_area_vals)), rep(median(MD_area_KBA_vals), length(MD_area_KBA_vals)),
          rep(median(LH_area_vals), length(LH_area_vals)), rep(median(LH_area_KBA_vals), length(LH_area_KBA_vals)),
          rep(median(VV_area_vals), length(VV_area_vals)), rep(median(VV_area_KBA_vals), length(VV_area_KBA_vals)),
          rep(median(BT_area_vals), length(BT_area_vals)), rep(median(BT_area_KBA_vals), length(BT_area_KBA_vals)),
          rep(median(HA_area_vals), length(HA_area_vals)), rep(median(HA_area_KBA_vals), length(HA_area_KBA_vals)))
val_df <- data.frame(nms, type, vals, meds)

val_plot <- ggplot(val_df, aes(x = vals, y = nms, fill = type)) +
  geom_violin(draw_quantiles = 0.5, 
              adjust = 0.2, #adjust changes the smoothness; lower is more faithful to the data
              scale = "width") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.text.y = element_text(face = "italic", size = 12),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  scale_fill_manual(values = c("#66B2FF", "#FFB266"),
                    labels = c("KBA mask", "species + area + weight"),
                    name = "Scenarios") +
  scale_y_discrete(limits = rev) +
  xlab("Site sensitivty")

#Kolmogorov-smirnoff tests - KBA vs no KBA
#a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution
#presence of ties a result of rounding
PM_KBA_area <- ks.test(PM_area_KBA_vals, PM_area_vals, alternative = "two.sided")
VG_KBA_area <- ks.test(VG_area_KBA_vals, VG_area_vals, alternative = "two.sided")
DG_KBA_area <- ks.test(DG_area_KBA_vals, DG_area_vals, alternative = "two.sided")
TB_KBA_area <- ks.test(TB_area_KBA_vals, TB_area_vals, alternative = "two.sided")
PL_KBA_area <- ks.test(PL_area_KBA_vals, PL_area_vals, alternative = "two.sided")
AM_KBA_area <- ks.test(AM_area_KBA_vals, AM_area_vals, alternative = "two.sided")
MF_KBA_area <- ks.test(MF_area_KBA_vals, MF_area_vals, alternative = "two.sided")
MD_KBA_area <- ks.test(MD_area_KBA_vals, MD_area_vals, alternative = "two.sided")
LH_KBA_area <- ks.test(LH_area_KBA_vals, LH_area_vals, alternative = "two.sided")
VV_KBA_area <- ks.test(VV_area_KBA_vals, VV_area_vals, alternative = "two.sided")
BT_KBA_area <- ks.test(BT_area_KBA_vals, BT_area_vals, alternative = "two.sided")
HA_KBA_area <- ks.test(HA_area_KBA_vals, HA_area_vals, alternative = "two.sided")

kbas_area <- tribble(
  ~species, ~statistic,
  "Pheidole megacephala", paste0(round(PM_KBA_area$statistic,4),"***"), 
  "Vespula germanica", paste0(round(VG_KBA_area$statistic,4),"***"), 
  "Digitonthophagus gazella", paste0(round(DG_KBA_area$statistic,4),"***"), 
  "Paratrechina longicornis", paste0(round(PL_KBA_area$statistic,4),"***"), 
  "Tetramorium bicarinatum", paste0(round(TB_KBA_area$statistic,4),"***"), 
  "Apis mellifera", paste0(round(AM_KBA_area$statistic,4),"***"), 
  "Monomorium floricola", paste0(round(MF_KBA_area$statistic,4),"***"), 
  "Monomorium destructor", paste0(round(MD_KBA_area$statistic,4),"***"), 
  "Linepithema humile", paste0(round(LH_KBA_area$statistic,4),"***"),
  "Vespula vulgaris", paste0(round(VV_KBA_area$statistic,4),"***"),
  "Bombus terrestris", paste0(round(BT_KBA_area$statistic,4),"***"), 
  "Heteronychus arator", paste0(round(HA_KBA_area$statistic,4),"***"), 
)

#Difference between IAS - species + area + weight
PM_VG_area <- ks.test(PM_area_vals, VG_area_vals, alternative = "two.sided")
PM_DG_area <- ks.test(PM_area_vals, DG_area_vals, alternative = "two.sided")
PM_TB_area <- ks.test(PM_area_vals, TB_area_vals, alternative = "two.sided")
PM_PL_area <- ks.test(PM_area_vals, PL_area_vals, alternative = "two.sided")
PM_AM_area <- ks.test(PM_area_vals, AM_area_vals, alternative = "two.sided")
PM_MF_area <- ks.test(PM_area_vals, MF_area_vals, alternative = "two.sided")
PM_MD_area <- ks.test(PM_area_vals, MD_area_vals, alternative = "two.sided")
PM_LH_area <- ks.test(PM_area_vals, LH_area_vals, alternative = "two.sided")
PM_VV_area <- ks.test(PM_area_vals, VV_area_vals, alternative = "two.sided")
PM_BT_area <- ks.test(PM_area_vals, BT_area_vals, alternative = "two.sided")
PM_HA_area <- ks.test(PM_area_vals, HA_area_vals, alternative = "two.sided")

VG_DG_area <- ks.test(VG_area_vals, DG_area_vals, alternative = "two.sided")
VG_TB_area <- ks.test(VG_area_vals, TB_area_vals, alternative = "two.sided")
VG_PL_area <- ks.test(VG_area_vals, PL_area_vals, alternative = "two.sided")
VG_AM_area <- ks.test(VG_area_vals, AM_area_vals, alternative = "two.sided")
VG_MF_area <- ks.test(VG_area_vals, MF_area_vals, alternative = "two.sided")
VG_MD_area <- ks.test(VG_area_vals, MD_area_vals, alternative = "two.sided")
VG_LH_area <- ks.test(VG_area_vals, LH_area_vals, alternative = "two.sided")
VG_VV_area <- ks.test(VG_area_vals, VV_area_vals, alternative = "two.sided")
VG_BT_area <- ks.test(VG_area_vals, BT_area_vals, alternative = "two.sided")
VG_HA_area <- ks.test(VG_area_vals, HA_area_vals, alternative = "two.sided")

DG_TB_area <- ks.test(DG_area_vals, TB_area_vals, alternative = "two.sided")
DG_PL_area <- ks.test(DG_area_vals, PL_area_vals, alternative = "two.sided")
DG_AM_area <- ks.test(DG_area_vals, AM_area_vals, alternative = "two.sided")
DG_MF_area <- ks.test(DG_area_vals, MF_area_vals, alternative = "two.sided")
DG_MD_area <- ks.test(DG_area_vals, MD_area_vals, alternative = "two.sided")
DG_LH_area <- ks.test(DG_area_vals, LH_area_vals, alternative = "two.sided")
DG_VV_area <- ks.test(DG_area_vals, VV_area_vals, alternative = "two.sided")
DG_BT_area <- ks.test(DG_area_vals, BT_area_vals, alternative = "two.sided")
DG_HA_area <- ks.test(DG_area_vals, HA_area_vals, alternative = "two.sided")

TB_PL_area <- ks.test(TB_area_vals, PL_area_vals, alternative = "two.sided")
TB_AM_area <- ks.test(TB_area_vals, AM_area_vals, alternative = "two.sided")
TB_MF_area <- ks.test(TB_area_vals, MF_area_vals, alternative = "two.sided")
TB_MD_area <- ks.test(TB_area_vals, MD_area_vals, alternative = "two.sided")
TB_LH_area <- ks.test(TB_area_vals, LH_area_vals, alternative = "two.sided")
TB_VV_area <- ks.test(TB_area_vals, VV_area_vals, alternative = "two.sided")
TB_BT_area <- ks.test(TB_area_vals, BT_area_vals, alternative = "two.sided")
TB_HA_area <- ks.test(TB_area_vals, HA_area_vals, alternative = "two.sided")

PL_AM_area <- ks.test(PL_area_vals, AM_area_vals, alternative = "two.sided")
PL_MF_area <- ks.test(PL_area_vals, MF_area_vals, alternative = "two.sided")
PL_MD_area <- ks.test(PL_area_vals, MD_area_vals, alternative = "two.sided")
PL_LH_area <- ks.test(PL_area_vals, LH_area_vals, alternative = "two.sided")
PL_VV_area <- ks.test(PL_area_vals, VV_area_vals, alternative = "two.sided")
PL_BT_area <- ks.test(PL_area_vals, BT_area_vals, alternative = "two.sided")
PL_HA_area <- ks.test(PL_area_vals, HA_area_vals, alternative = "two.sided")

AM_MF_area <- ks.test(AM_area_vals, MF_area_vals, alternative = "two.sided")
AM_MD_area <- ks.test(AM_area_vals, MD_area_vals, alternative = "two.sided")
AM_LH_area <- ks.test(AM_area_vals, LH_area_vals, alternative = "two.sided")
AM_VV_area <- ks.test(AM_area_vals, VV_area_vals, alternative = "two.sided")
AM_BT_area <- ks.test(AM_area_vals, BT_area_vals, alternative = "two.sided")
AM_HA_area <- ks.test(AM_area_vals, HA_area_vals, alternative = "two.sided")

MF_MD_area <- ks.test(MF_area_vals, MD_area_vals, alternative = "two.sided")
MF_LH_area <- ks.test(MF_area_vals, LH_area_vals, alternative = "two.sided")
MF_VV_area <- ks.test(MF_area_vals, VV_area_vals, alternative = "two.sided")
MF_BT_area <- ks.test(MF_area_vals, BT_area_vals, alternative = "two.sided")
MF_HA_area <- ks.test(MF_area_vals, HA_area_vals, alternative = "two.sided")

MD_LH_area <- ks.test(MD_area_vals, LH_area_vals, alternative = "two.sided")
MD_VV_area <- ks.test(MD_area_vals, VV_area_vals, alternative = "two.sided")
MD_BT_area <- ks.test(MD_area_vals, BT_area_vals, alternative = "two.sided")
MD_HA_area <- ks.test(MD_area_vals, HA_area_vals, alternative = "two.sided")

LH_VV_area <- ks.test(LH_area_vals, VV_area_vals, alternative = "two.sided")
LH_BT_area <- ks.test(LH_area_vals, BT_area_vals, alternative = "two.sided")
LH_HA_area <- ks.test(LH_area_vals, HA_area_vals, alternative = "two.sided")

VV_BT_area <- ks.test(VV_area_vals, BT_area_vals, alternative = "two.sided")
VV_HA_area <- ks.test(VV_area_vals, HA_area_vals, alternative = "two.sided")

BT_HA_area <- ks.test(BT_area_vals, HA_area_vals, alternative = "two.sided")

v <- c(0,PM_VG_area$statistic,PM_DG_area$statistic,PM_TB_area$statistic,PM_PL_area$statistic,PM_AM_area$statistic,PM_MF_area$statistic,PM_MD_area$statistic,PM_LH_area$statistic,PM_VV_area$statistic,PM_BT_area$statistic,PM_HA_area$statistic,
       0,0,VG_DG_area$statistic,VG_TB_area$statistic,VG_PL_area$statistic,VG_AM_area$statistic,VG_MF_area$statistic,VG_MD_area$statistic,VG_LH_area$statistic,VG_VV_area$statistic,VG_BT_area$statistic,VG_HA_area$statistic,
       0,0,0,DG_TB_area$statistic,DG_PL_area$statistic,DG_AM_area$statistic,DG_MF_area$statistic,DG_MD_area$statistic,DG_LH_area$statistic,DG_VV_area$statistic,DG_BT_area$statistic,DG_HA_area$statistic,
       0,0,0,0,TB_PL_area$statistic,TB_AM_area$statistic,TB_MF_area$statistic,TB_MD_area$statistic,TB_LH_area$statistic,TB_VV_area$statistic,TB_BT_area$statistic,TB_HA_area$statistic,
       0,0,0,0,0,PL_AM_area$statistic,PL_MF_area$statistic,PL_MD_area$statistic,PL_LH_area$statistic,PL_VV_area$statistic,PL_BT_area$statistic,PL_HA_area$statistic,
       0,0,0,0,0,0,AM_MF_area$statistic,AM_MD_area$statistic,AM_LH_area$statistic,AM_VV_area$statistic,AM_BT_area$statistic,AM_HA_area$statistic,
       0,0,0,0,0,0,0,MF_MD_area$statistic,MF_LH_area$statistic,MF_VV_area$statistic,MF_BT_area$statistic,MF_HA_area$statistic,
       0,0,0,0,0,0,0,0,MD_LH_area$statistic,MD_VV_area$statistic,MD_BT_area$statistic,MD_HA_area$statistic,
       0,0,0,0,0,0,0,0,0,LH_VV_area$statistic,LH_BT_area$statistic,LH_HA_area$statistic,
       0,0,0,0,0,0,0,0,0,0,VV_BT_area$statistic,VV_HA_area$statistic,
       0,0,0,0,0,0,0,0,0,0,0,BT_HA_area$statistic,
       0,0,0,0,0,0,0,0,0,0,0,0)

tm <- matrix(v, nrow = 12, ncol = 12)
rownames(tm) <- c(":italic(Pheidole~~megacephala)", 
                  ":italic(Vespula~~germanica)", 
                  ":italic(Digitonthophagus~~gazella)", 
                  ":italic(Tetramorium~~bicarinatum)", 
                  ":italic(Paratrechina~~longicornis)",
                  ":italic(Apis~~mellifera)",
                  ":italic(Monomorium~~floricola)",
                  ":italic(Monomorium~~destructor)",
                  ":italic(Linepithema~~humile)",
                  ":italic(Vespula~~vulgaris)",
                  ":italic(Bombus~~terrestris)",
                  ":italic(Heteronychus~~arator)")
colnames(tm) <- c("Pm", "Vg", "Dg", "Tb", "Pl", "Am", "Mf", "Md", "Lh", "Vv", "Bt", "Ha")
corrplot::corrplot(tm, type = "lower", method = "color", 
                   cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
                   tl.srt = 0, tl.col = "black", tl.cex = 0.8, addCoef.col = "black",
                   mar = c(0,0,0,0))

# #Alternative
# rv <- rev(v)
# tmr <- matrix(rv, nrow = 12, ncol = 12)
# rownames(tmr) <- rev(c(":italic(Pheidole~~megacephala)", 
#                   ":italic(Vespula~~germanica)", 
#                   ":italic(Digitonthophagus~~gazella)", 
#                   ":italic(Tetramorium~~bicarinatum)", 
#                   ":italic(Paratrechina~~longicornis)",
#                   ":italic(Apis~~mellifera)",
#                   ":italic(Monomorium~~floricola)",
#                   ":italic(Monomorium~~destructor)",
#                   ":italic(Linepithema~~humile)",
#                   ":italic(Vespula~~vulgaris)",
#                   ":italic(Bombus~~terrestris)",
#                   ":italic(Heteronychus~~arator)"))
# colnames(tm) <- rev(c("Pm", "Vg", "Dg", "Tb", "Pl", "Am", "Mf", "Md", "Lh", "Vv", "Bt", "Ha"))
# corrplot::corrplot(tm, type = "upper", method = "color", 
#                    cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
#                    tl.srt = 0, tl.col = "black", tl.cex = 0.8, addCoef.col = "black",
#                    mar = c(0,0,0,0))

#Proportion difference (KBA vs no KBA) in number of top sensitive sites
#Top two (i.e. >= 0.98 sensitivity)
props <- c(1, 0.98, 0.95, 0.90, 0.75, 0.50, 0.25, 0.00)
diffs <- c()

PM_area_prop <- multi_props(PM_area_vals, props)
PM_area_KBA_prop <- multi_props(PM_area_KBA_vals, props)
PM_area_RAN_prop <- multi_props(PM_area_RAN_vals, props)

VG_area_prop <- multi_props(VG_area_vals, props)
VG_area_KBA_prop <- multi_props(VG_area_KBA_vals, props)
VG_area_RAN_prop <- multi_props(VG_area_RAN_vals, props)

DG_area_prop <- multi_props(DG_area_vals, props)
DG_area_KBA_prop <- multi_props(DG_area_KBA_vals, props)
DG_area_RAN_prop <- multi_props(DG_area_RAN_vals, props)

TB_area_prop <- multi_props(TB_area_vals, props)
TB_area_KBA_prop <- multi_props(TB_area_KBA_vals, props)
TB_area_RAN_prop <- multi_props(TB_area_RAN_vals, props)

PL_area_prop <- multi_props(PL_area_vals, props)
PL_area_KBA_prop <- multi_props(PL_area_KBA_vals, props)
PL_area_RAN_prop <- multi_props(PL_area_RAN_vals, props)

AM_area_prop <- multi_props(AM_area_vals, props)
AM_area_KBA_prop <- multi_props(AM_area_KBA_vals, props)
AM_area_RAN_prop <- multi_props(AM_area_RAN_vals, props)

MF_area_prop <- multi_props(MF_area_vals, props)
MF_area_KBA_prop <- multi_props(MF_area_KBA_vals, props)
MF_area_RAN_prop <- multi_props(MF_area_RAN_vals, props)

MD_area_prop <- multi_props(MD_area_vals, props)
MD_area_KBA_prop <- multi_props(MD_area_KBA_vals, props)
MD_area_RAN_prop <- multi_props(MD_area_RAN_vals, props)

LH_area_prop <- multi_props(LH_area_vals, props)
LH_area_KBA_prop <- multi_props(LH_area_KBA_vals, props)
LH_area_RAN_prop <- multi_props(LH_area_RAN_vals, props)

VV_area_prop <- multi_props(VV_area_vals, props)
VV_area_KBA_prop <- multi_props(VV_area_KBA_vals, props)
VV_area_RAN_prop <- multi_props(VV_area_RAN_vals, props)

BT_area_prop <- multi_props(BT_area_vals, props)
BT_area_KBA_prop <- multi_props(BT_area_KBA_vals, props)
BT_area_RAN_prop <- multi_props(BT_area_RAN_vals, props)

HA_area_prop <- multi_props(HA_area_vals, props)
HA_area_KBA_prop <- multi_props(HA_area_KBA_vals, props)
HA_area_RAN_prop <- multi_props(HA_area_RAN_vals, props)

Total <- as.data.frame(rbind(PM_area_prop,VG_area_prop,DG_area_prop,TB_area_prop,PL_area_prop,AM_area_prop,MF_area_prop,MD_area_prop,LH_area_prop,VV_area_prop,BT_area_prop,HA_area_prop,
                             PM_area_KBA_prop, VG_area_KBA_prop, DG_area_KBA_prop, TB_area_KBA_prop, PL_area_KBA_prop,AM_area_KBA_prop,MF_area_KBA_prop,MD_area_KBA_prop,LH_area_KBA_prop,VV_area_KBA_prop,BT_area_KBA_prop,HA_area_KBA_prop,
                             PM_area_RAN_prop, VG_area_RAN_prop, DG_area_RAN_prop, TB_area_RAN_prop, PL_area_RAN_prop,AM_area_RAN_prop,MF_area_RAN_prop,MD_area_RAN_prop,LH_area_RAN_prop,VV_area_RAN_prop,BT_area_RAN_prop,HA_area_RAN_prop))
Type <- c(rep("Unmasked", 12), rep("Masked", 12), rep("Random", 12))
Species <- rep(c("P. megacephala", "V. germanica", "D. gazella",
                 "T. bicarinatum", "P. longicornis", "A. mellifera",
                 "M. floricola", "M. destructor", "L. humile",
                 "V. vulgaris", "B. terrestris", "H. arator"), 3)
Total <- cbind(Species, Type,Total)
Total <- Total %>%
  as_tibble() %>%
  mutate(Species = factor(Species)) %>%
  mutate(Type = factor(Type)) %>%
  pivot_longer(cols = !Species & !Type, 
               names_to = "SiteSensitivity", 
               values_to = "DistributionCoverage") %>%
  mutate(SiteSensitivity = as.double(SiteSensitivity))
Total_min <- Total %>%
  dplyr::filter(SiteSensitivity == 1 | SiteSensitivity == 0.98)
df <- data.frame(xmin = c(1,0.98,0.95,0.9,0.75,0.5,0.25),
                 xmax = c(0.98,0.95,0.9,0.75,0.5,0.25,0),
                 ymin = rep(-Inf,7),
                 ymax = rep(Inf,7))
p <- ggplot(Total, 
            aes(x = SiteSensitivity, 
                y = DistributionCoverage, 
                colour = Species,
                group = interaction(Species,Type))) +
  geom_line(stat = "identity", aes(linetype = Type), size = 0.8) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) +
  ylab("Distribution Coverage (%)") +
  xlab("Site Sensitivity") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0))
p2 <-  p + 
  geom_rect(data = df, aes(xmin = rev(xmin), xmax = rev(xmax), ymin=ymin, ymax=ymax),
            fill = leg$colors, alpha = 0.2, inherit.aes = F) +
  font("legend.text", face = "italic")

#For inset
ins <- ggplot(Total_min, 
              aes(x = SiteSensitivity, 
                  y = DistributionCoverage, 
                  colour = Species,
                  group = interaction(Species,Type))) +
  geom_line(stat = "identity", aes(linetype = Type),size = 0.8) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  #scale_size_manual(values = c(1,1,1)) +
  scale_colour_manual(values = pnw_palette("Starfish", 12)) +
  ylab("Distribution Coverage (%)") +
  xlab("Site Sensitivity") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = alpha(leg$colors[7],0.2)),
        axis.text = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank()) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0.98,1)) +
  scale_y_continuous(expand = c(0,0))

#An accompanying figure could be showing numerical difference between non-kba and kba
#Take differences, e.g. Here negative values means higher distribution coverage for KBAs
pm_area_diff <- PM_area_prop - PM_area_KBA_prop
vg_area_diff <- VG_area_prop - VG_area_KBA_prop
dg_area_diff <- DG_area_prop - DG_area_KBA_prop
tb_area_diff <- TB_area_prop - TB_area_KBA_prop
pl_area_diff <- PL_area_prop - PL_area_KBA_prop
am_area_diff <- AM_area_prop - AM_area_KBA_prop
mf_area_diff <- MF_area_prop - MF_area_KBA_prop
md_area_diff <- MD_area_prop - MD_area_KBA_prop
lh_area_diff <- LH_area_prop - LH_area_KBA_prop
vv_area_diff <- VV_area_prop - VV_area_KBA_prop
bt_area_diff <- BT_area_prop - BT_area_KBA_prop
ha_area_diff <- HA_area_prop - HA_area_KBA_prop

#Then just follow through like before but using these values
Total_diff_area <- data.frame(rbind(pm_area_diff, vg_area_diff, dg_area_diff, tb_area_diff, 
                               pl_area_diff, am_area_diff, mf_area_diff, md_area_diff, 
                               lh_area_diff, vv_area_diff, bt_area_diff, ha_area_diff))
Species <- c("Pheidole megacephala", "Vespula germanica", "Digitonthophagus gazella",
             "Tetramorium bicarinatum", "Paratrechina longicornis",
             "Apis mellifera", "Monomorium floricola", "Monomorium destructor",
             "Linepithema humile", "Vespula vulgaris", "Bombus terrestris", 
             "Heteronychus arator")
Total_diff_area <- data.frame(Species, Total_diff_area)
names(Total_diff_area) <- c("Species", "1.00","0.98","0.95","0.90","0.75","0.50","0.25", "0.00")
Total_diff_area <- Total_diff_area %>%
  as_tibble() %>%
  mutate(Species = factor(Species)) %>%
  pivot_longer(cols = !Species, 
               names_to = "SiteSensitivity", 
               values_to = "ScenarioDifference") %>%
  mutate(SiteSensitivity = as.double(SiteSensitivity))

Total_diff_area <- Total_diff_area %>% 
  group_by(SiteSensitivity) %>% 
  mutate(mu = mean(ScenarioDifference)) %>% 
  mutate(sd = sqrt(var(ScenarioDifference)))

ggplot(Total_diff_area, 
       aes(x = SiteSensitivity, 
           y = ScenarioDifference, 
           colour = Species)) +
  geom_line(stat = "identity", size = 0.8) +
  scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  geom_hline(yintercept = 0) +
  geom_rect(data = df, aes(xmin = rev(xmin), xmax = rev(xmax), ymin=ymin, ymax=ymax),
            fill = leg$colors, alpha = 0.2, inherit.aes = F) +
  font("legend.text", face = "italic") +
  ylab("Scenario difference in distribution coverage") +
  xlab("Site sensitivity")

#Mean difference among species
ggplot(Total_diff_area, 
       aes(x = SiteSensitivity, 
           y = ScenarioDifference)) +
  geom_smooth(col = "black", fill = "black", alpha = 0.1) +
  geom_smooth(data = Total_diff, aes(x = SiteSensitivity, 
                                     y = ScenarioDifference), 
              col = "purple", fill = "purple", alpha = 0.1) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  geom_hline(yintercept = 0) +
  geom_rect(data = df, aes(xmin = rev(xmin), xmax = rev(xmax), ymin=ymin, ymax=ymax),
            fill = leg$colors, alpha = 0.2, inherit.aes = F) +
  font("legend.text", face = "italic") +
  ylab("Scenario difference in distribution coverage") +
  xlab("Site sensitivity")


#Multi-IAS priority sites----

CAZ_wgt_var_ras <- rast("C:/Users/david/Documents/PhD/Chapter_3/SpatialData/CAZ_wgt_var_ras.tif")
ias_sum_sp <- rast("C:/Users/david/Documents/PhD/Chapter_3/SpatialData/ias_richness_sp.tif")

#Need an empty new raster
ext <- extent(ias_sum_sp_ar)
crs <- "+proj=longlat +datum=WGS84 +no_defs"
res <- 0.04166667
vals <- 0
new_raster_0.5 <- raster(ext = ext, crs = crs, resolution = res, vals = vals)
new_raster_0.5 <- mask(new_raster_0.5, ias_sum_sp_ar)

new_raster_0.7 <- raster(ext = ext, crs = crs, resolution = res, vals = vals)
new_raster_0.7 <- mask(new_raster_0.7, ias_sum_sp_ar)

#Weighted species only sensitive sites (maybe do species and area; focus on this for paper)
#Arbitrary levels of IAS richness and site sensitivity (mid levels of both)
ias_strat_1 <- which(ias_sum_sp_ar[] < 6)
caz_strat_1 <- which(CAZ_area_wgt_var_ras[] < 0.5)
strat_1 <- intersect(ias_strat_1, caz_strat_1)

ias_strat_2 <- which(ias_sum_sp_ar[] >= 6)
caz_strat_2 <- which(CAZ_area_wgt_var_ras[] < 0.5)
strat_2 <- intersect(ias_strat_2, caz_strat_2)

ias_strat_3 <- which(ias_sum_sp_ar[] < 6)
caz_strat_3 <- which(CAZ_area_wgt_var_ras[] >= 0.5)
strat_3 <- intersect(ias_strat_3, caz_strat_3)

ias_strat_4 <- which(ias_sum_sp_ar[] >= 6)
caz_strat_4 <- which(CAZ_area_wgt_var_ras[] >= 0.5)
strat_4 <- intersect(ias_strat_4, caz_strat_4)

#Management strategy map
new_raster_0.5[strat_1] <- 1
new_raster_0.5[strat_2] <- 2
new_raster_0.5[strat_3] <- 3
new_raster_0.5[strat_4] <- 4
new_raster_0.5[new_raster_0.5 == 0] <- NA
names(new_raster_0.5) <- "strategy"

ias_strat_1 <- which(ias_sum_sp_ar[] < 3)
caz_strat_1 <- which(CAZ_area_wgt_var_ras[] < 0.7)
strat_1 <- intersect(ias_strat_1, caz_strat_1)

ias_strat_2 <- which(ias_sum_sp_ar[] >= 3)
caz_strat_2 <- which(CAZ_area_wgt_var_ras[] < 0.7)
strat_2 <- intersect(ias_strat_2, caz_strat_2)

ias_strat_3 <- which(ias_sum_sp_ar[] < 3)
caz_strat_3 <- which(CAZ_area_wgt_var_ras[] >= 0.7)
strat_3 <- intersect(ias_strat_3, caz_strat_3)

ias_strat_4 <- which(ias_sum_sp_ar[] >= 3)
caz_strat_4 <- which(CAZ_area_wgt_var_ras[] >= 0.7)
strat_4 <- intersect(ias_strat_4, caz_strat_4)

#Management strategy map
new_raster_0.7[strat_1] <- 1
new_raster_0.7[strat_2] <- 2
new_raster_0.7[strat_3] <- 3
new_raster_0.7[strat_4] <- 4
new_raster_0.7[new_raster_0.7 == 0] <- NA
names(new_raster_0.7) <- "strategy"

rm(ias_strat_1, ias_strat_2, ias_strat_3, ias_strat_4,
   caz_strat_1, caz_strat_2, caz_strat_3, caz_strat_4,
   strat_1, strat_2, strat_3, strat_4)

#need row number
daf <- data.frame(susc = values(ias_sum_sp_ar), 
                  sens = values(CAZ_area_wgt_var_ras), 
                  strategy = values(new_raster_0.5),
                  cellID = seq(1:length(values(ias_sum_sp_ar))))

sd1 <- get_shade(daf, "susc", "sens", 1, c(0,1), c(0,1)) 
sd2 <- get_shade(daf, "susc", "sens", 2, c(0,1), c(1,0))
sd3 <- get_shade(daf, "susc", "sens", 3, c(1,0), c(0,1))
sd4 <- get_shade(daf, "susc", "sens", 4, c(1,0), c(1,0))

sd <- sd1 %>%
  bind_rows(sd2) %>%
  bind_rows(sd3) %>%
  bind_rows(sd4) %>%
  arrange(cellID)


daf_sd <- daf %>% left_join(sd, by = c("cellID" = "cellID", 
                                       "susc" = "susc", 
                                       "sens" = "sens",
                                       "strategy" = "strategy"))

#Area in QLD
xmn <- dm2dd(148, 27, 10.58)
xmx <- dm2dd(151, 23, 05.29)
ymn <- dm2dd(23, 40, 34.46, "S")
ymx <- dm2dd(20, 03, 21.49, "S")

#Extent object
qld_ex <- extent(c(xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx))
qld_cells <- cellsFromExtent(new_raster_0.5, qld_ex)

ras_st <- st_as_stars(new_raster_0.5)
ras_sf <- st_as_sf(ras_st)
daf_ras <- ras_sf %>%
  bind_cols(sd) %>%
  dplyr::select(-"strategy...5") %>%
  rename(strategy = "strategy...1") %>%
  mutate(susc_res = scales::rescale(susc, to = c(0,1)))

fill_weight_ras_sf <- daf_ras %>%
  dplyr::select(weight) %>%
  st_drop_geometry()

fill_strat_ras_sf <- daf_ras %>%
  dplyr::select(strategy) %>%
  st_drop_geometry()



strat_map <- ggplot()+
  geom_sf(data = ras_sf, 
          aes(fill = factor(fill_strat_ras_sf[,1]),
              alpha = fill_weight_ras_sf[,1]), 
          color = NA, 
          show.legend = T) + 
  scale_fill_brewer(palette = "Dark2", 
                    type = "qual",
                    name = "Strategy") +
  # scale_alpha_continuous(range = c(0,1),
  #                        name = "Conviction") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  geom_rect(aes(xmin = xmn, 
                xmax = xmx, 
                ymin = ymn, 
                ymax = ymx), 
            colour = "black", 
            fill = NA,
            size = 2)

ggsave(filename = "strat_map_spec_area.pdf", 
       plot = strat_map, 
       device = "pdf", 
       width = 5840,
       height = 4132,
       units = "px",
       path = file.path("Outcome", "pdf"))

#Filter QLD cells
qld_daf <- daf_ras %>%
  filter(cellID %in% qld_cells)

fill_weight_ras_sf_qld <- qld_daf %>%
  dplyr::select(weight) %>%
  st_drop_geometry()

fill_strat_ras_sf_qld <- qld_daf %>%
  dplyr::select(strategy) %>%
  st_drop_geometry()

strat_map_qld <- ggplot()+
  geom_sf(data = qld_daf, 
          aes(fill = factor(fill_strat_ras_sf_qld[,1]),
              alpha = fill_weight_ras_sf_qld[,1]), 
          color = NA, 
          show.legend = T) + 
  scale_fill_brewer(palette = "Dark2", 
                    type = "qual",
                    name = "Strategy",
                    labels = c("Maintain", "Containment", "Prevention", "Eradication")) +
  scale_alpha_continuous(name = "Conviction") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))

gg_inset_map1 <- ggdraw() +
  draw_plot(strat_map_qld) +
  draw_plot(strat_map, x = 0.45, y = 0.65, width = 0.3, height = 0.3)

ggsave(filename = "strat_map_qld.pdf", 
       plot = gg_inset_map1, 
       device = "pdf", 
       width = 5840,
       height = 4132,
       units = "px",
       path = file.path("Outcome", "pdf"))


#rich_sp <- overlay(ias_sum_sp, CAZ_wgt_var_ras, fun = function(x,y){ias_wght(x,y)})

#Can I get a weighting on the -1 to 1 scale
#(1,1) top right, (-1,1) top left, (-1,-1) bottom left, (1,-1) bottom right

# df <- data.frame(susc = values(ias_sum_sp),
#                  sens = values(CAZ_wgt_var_ras))
# df <- df %>%
#   mutate(sb_weight = ias_wght(susc,sens)) %>%
#   mutate(susc_sc = scales::rescale(susc, to = c(0,1))) %>%
#   mutate(strategy = case_when(susc_sc < 0.5454545 & sens < 0.5 ~ "Strat1",
#                               susc_sc >= 0.5454545 & sens < 0.5 ~ "Strat2",
#                               susc_sc < 0.5454545 & sens >= 0.5 ~ "Strat3",
#                               susc_sc >= 0.5454545 & sens >= 0.5 ~ "Strat4")) %>%
#   mutate(strategy = as.fator(strategy))
#   #mutate(conf = case_when(susc_sc <= 0.1 & sens <= 0.1 ~ ))
# 
# #could filter by strategy, put all values 0-1 scale, find elegant pairing
# #this could be to create a shading variable (elegant means unordered which is fine for this purpose)
# new_df <- df %>%
#   group_by(strategy) %>%
#   drop_na(strategy) %>%
#   slice_sample(n = 500)


# #Strategy 1
# strat_1_df <- get_shade(new_df, "susc", "sens", "Strat1")
# 
# #Strategy 2
# strat_2_df <- get_shade(new_df, "susc", "sens", "Strat2")
# 
# #Strategy 3
# strat_3_df <- get_shade(new_df, "susc", "sens", "Strat3")
# 
# #Strategy 4
# strat_4_df <- get_shade(new_df, "susc", "sens", "Strat4")
# 
# #Combine
# new_df_wgt <- strat_1_df %>%
#   bind_rows(strat_2_df) %>%
#   bind_rows(strat_3_df) %>%
#   bind_rows(strat_4_df)

dat25 <- circleFun(c(0.5,0.5),0.25,npoints = 100)
dat50 <- circleFun(c(0.5,0.5),0.50,npoints = 100)
dat75 <- circleFun(c(0.5,0.5),0.75,npoints = 100)
dat100 <- circleFun(c(0.5,0.5),1,npoints = 100)
dat120 <- circleFun(c(0.5,0.5),2,npoints = 100)

bg_cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))
#Antarctica = c("#1C2336", "#3B475D", "#5A6C85", "#7A91AD", "#97AEC8", "#B1C2D4","#CCD5E1", "#E7EAEE") #Just install my package
#bg_cols <- colorRampPalette(rev(Antarctica))

ll <- make_gradient(
  deg = 135, n = 500, cols = bg_cols(20)
)

lh <- make_gradient(
  deg = 225, n = 500, cols = bg_cols(20)
)

hh <- make_gradient(
  deg = 315, n = 500, cols = bg_cols(20)
)

hl <- make_gradient(
  deg = 45, n = 500, cols = bg_cols(20)
)

ggplot(qld_daf, aes(x = sens, y = susc_res)) +
  annotation_custom(
    grob = ll, xmin = 0.0, xmax = 0.5, ymin = 0.0, ymax = 0.5
  ) +
  annotation_custom(
    grob = lh, xmin = 0.0, xmax = 0.5, ymin = 0.5, ymax = 1.0
  ) +
  annotation_custom(
    grob = hh, xmin = 0.5, xmax = 1.0, ymin = 0.5, ymax = 1.0
  ) +
  annotation_custom(
    grob = hl, xmin = 0.5, xmax = 1.0, ymin = 0.0, ymax = 0.5
  ) +
  geom_point(aes(fill = factor(strategy)), shape = 21, size = 4) +
  #scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(breaks = c(1,2,3,4),
                    values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
                    labels = c("Maintain", "Containment", "Prevention", "Eradication"),
                    name = "Strategy") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_y_continuous(labels = c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7",
                                "0.8","0.9","1.0"),
                     breaks = c(seq(0,1, by = 0.1)),
                     limits = c(0,1)) +
  ylab("Species richness (scaled)") +
  xlab("Site sensitivity") +
  geom_hline(yintercept = 0.5, col = "black") +
  geom_vline(xintercept = 0.5, col = "black") +
  geom_path(data = dat25, mapping = aes(x,y)) +
  geom_path(data = dat50, mapping = aes(x,y)) +
  geom_path(data = dat75, mapping = aes(x,y)) +
  geom_path(data = dat100, mapping = aes(x,y))
