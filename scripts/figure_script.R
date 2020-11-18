# Zinc Manuscript Figures

###############################
# IMPORTANT!!!!!
# Set working directory
setwd("~/Dropbox/AndersenLab/LabFolders/Katie/git/zinc_manuscript/")
###############################

options(stringsAsFactors=FALSE)

library(tidyverse)
library(linkagemapping)
library(cegwas2)
library(ggtree)
library(phangorn)
library(ape)
# library(grid)
# library(gridExtra)
library(ggplotify)
library(easysorter)

tsize <- 12


##########################################
#           FUNCTIONS                    #
##########################################

load("data/extra_data/gene_annotations.Rda")
load("data/extra_data/FileS9_eqtlmap.Rda")
load("data/extra_data/eqtl_probes.Rda")
load("data/extra_data/FileS8_eqtlpheno.Rda")

# linkage mapping function 
plot_lods <- function(map, tsize = 12) {
    map1 <- map %>% 
        dplyr::group_by(marker, condtrt) %>% 
        dplyr::filter(lod == max(lod))
    cis <- map %>% 
        dplyr::group_by(marker, condtrt) %>% 
        dplyr::mutate(maxlod = max(lod)) %>%
        dplyr::group_by(iteration, condtrt) %>% 
        dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n = 1))
    
    totalmap <- NULL
    if(nrow(cis) > 0) {
        for(i in unique(cis$condtrt)) {
            drugci <- cis %>%
                dplyr::filter(condtrt == i)
            drugmap <- map1 %>%
                dplyr::filter(condtrt == i)
            map2 <- linkagemapping:::cidefiner(drugci, drugmap)
            totalmap <- rbind(totalmap, map2)
        }
        
        totalmap$condtrt <- gsub("_", "\n", totalmap$condtrt)
        cis$condtrt <- gsub("_", "\n", cis$condtrt)
        
        ggplot2::ggplot(totalmap) + 
            ggplot2::aes(x = pos/1e+06,y = lod) + 
            ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = "skyblue", alpha = 0.5) +
            ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                                size = tsize/7, show.legend = FALSE, color = "red") + 
            ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.2 * maxlod), label = paste0(100 *round(var_exp, digits = 3), "%")),
                               size = tsize/5, colour = "black", hjust = "inward") +
            # ggrepel::geom_text_repel(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.5 * maxlod), label = paste0(100 *round(var_exp, digits = 3), "%")), 
            #                    size = tsize/5, colour = "black") + 
            ggplot2::geom_line(size = tsize/25, alpha = 0.85) +
            ggplot2::facet_grid(~ chr, scales = "free", space = "free_x") +
            ggplot2::labs(x = "Genomic position (Mb)", y = "LOD") +
            ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
            # ggplot2::scale_x_continuous(expand = c(0,0)) + 
            # ggplot2::scale_y_continuous(expand = expand_scale(mult=c(0,0.1))) +
            ggplot2::theme_bw(tsize) +
            ggplot2::theme(
                axis.text = ggplot2::element_text(color = "black", face = "bold"),
                axis.title = ggplot2::element_text(face = "bold", color = "black"),
                strip.text = ggplot2::element_text(face = "bold", color = "black"),
                plot.title = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(color = NA, size = 0.6))
    } else {
        totalmap <- map1
        totalmap$condtrt <- gsub("_", "\n", totalmap$condtrt)
        
        ggplot2::ggplot(totalmap) + 
            ggplot2::aes(x = pos/1e+06,y = lod) + 
            ggplot2::geom_line(size = tsize/25, alpha = 0.85) +
            ggplot2::facet_grid( ~ chr, scales = "free", space = "free_x") +
            ggplot2::labs(x = "Genomic position (Mb)", y = "LOD") +
            ggplot2::theme_bw(tsize) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(color = "black", face = "bold"),
                axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
                axis.title.x = ggplot2::element_text(face = "bold", color = "black"),
                axis.title.y = ggplot2::element_text(face = "bold", color = "black"),
                strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
                strip.text.y = ggplot2::element_text(face = "bold", color = "black"),
                plot.title = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(color = NA, size = 0.6))
    }
}

# function for PxG 
plot_pxg <- function(cross, map, tsize = 12, yaxis = "Relative animal length") {
    
    # get unique QTL peaks
    peaks <- map %>% 
        dplyr::group_by(iteration, condition) %>% 
        dplyr::filter(!is.na(var_exp)) %>% 
        dplyr::do(head(., n = 1))
    
    # clean the markers and column names
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    colnames(cross$pheno) <- stringr::str_replace(colnames(cross$pheno), "\\.", "_")
    
    # get only the traits of interest
    pheno <- cross$pheno %>% 
        dplyr::select(dplyr::one_of(map$condtrt))
    
    # get the genotype for the RIAILs and add to pheno
    geno <- data.frame(linkagemapping:::extract_genotype(cross)) %>% 
        dplyr::select(which(colnames(.) %in% uniquemarkers)) %>% 
        data.frame(., pheno)
    
    # reorder data and change -1 and 1 to N2 and CB and plot!
    df <- geno %>%
        tidyr::gather(marker, genotype, -dplyr::one_of(map$condtrt)) %>%
        dplyr::mutate(genotype = dplyr::case_when(genotype == -1 ~ "N2",
                                                  genotype == 1 ~ "CB4856",
                                                  TRUE ~ "NA")) %>%
        tidyr::gather(trait, phenotype, dplyr::one_of(map$condtrt)) %>%
        dplyr::mutate(genotype = factor(genotype, levels = c("N2", "CB4856"), labels= c("N2", "CB4856"))) %>%
        dplyr::left_join(peaks, by = c("marker", "trait" = "condtrt")) %>%
        tidyr::drop_na(lod) %>%
        dplyr::mutate(marker = stringr::str_replace(marker, "_", ":"),
                      trait = stringr::str_split_fixed(trait, "_", 2)[,2],
                      pos = as.numeric(stringr::str_split_fixed(marker, ":", 2)[,2])) %>%
        dplyr::arrange(chr, pos) %>%
        tidyr::drop_na(genotype)
    df %>%
        ggplot2::ggplot(.) + 
        ggplot2::aes(x = genotype, y = phenotype) +
        ggplot2::geom_jitter(width = 0.1, size = 0.07, alpha = 0.5) + 
        ggplot2::geom_boxplot(ggplot2::aes(fill = genotype, alpha = 0.5), size = 0.2, outlier.shape = NA) + 
        ggplot2::scale_fill_manual(values = c(`N2` = "orange", `CB4856` = "blue")) + 
        ggplot2::facet_grid(~factor(marker, unique(df$marker)), scales = "free") + 
        ggplot2::theme_bw(tsize) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", color = "black"), 
                       axis.text.y = ggplot2::element_text(face = "bold", color = "black"), 
                       axis.title.x = ggplot2::element_text(face = "bold", color = "black", vjust = -0.3),
                       axis.title.y = ggplot2::element_text(face = "bold", color = "black"), 
                       strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
                       strip.text.y = ggplot2::element_text(face = "bold", color = "black"),
                       # axis.title.x = element_blank(),
                       plot.title = ggplot2::element_blank(), 
                       legend.position = "none", 
                       panel.grid = element_blank(),
                       panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
        ggplot2::labs(x = "Genotype at QTL", y = yaxis)
}

# plot NIL pheno and geno with stats
plot_nil <- function(phenodf, genodf, statdf, strains, chr, tsize = 12, ylab = "Animal length") {
    
    # plot phenotype
    pheno <- phenodf %>%
        dplyr::group_by(strain, condition) %>%
        dplyr::mutate(phen = max(phenotype) + 0.1) %>%
        dplyr::ungroup() %>%
        dplyr::full_join(statdf, by = c("strain", "condition", "trait")) %>%
        dplyr::mutate(strain = factor(strain, 
                                      levels = rev(strains))) %>%
        dplyr::filter(!is.na(strain),
                      strain %in% strains) %>%
        ggplot(.) +
        aes(x = strain, y = phenotype, fill = strain_fill) +
        geom_jitter(width = 0.1, size = 0.05) +
        geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.2) +
        ggplot2::geom_text(aes(label = sig, y = phen, color = groups), size = tsize/4, angle = -90) +
        scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
        scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
        theme_bw(tsize) +
        theme(axis.text.x = element_text(face="bold", color="black"),
              axis.title.x = element_text(face="bold", color="black"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text = element_text(face = "bold", color = "black"),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) +
        coord_flip() +
        labs(x = " ", y = ylab)  +
        facet_grid(~condition)
    
    # plot chrV genotype
    chrgeno <- genodf %>%
        dplyr::filter(chrom == chr,
                      sample %in% strainset) %>%
        dplyr::mutate(chrom = paste0("chr", chr),
                      sample = factor(sample, levels = strainset)) %>%
        dplyr::distinct() %>%
        ggplot(.)+
        geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
        facet_grid(~chrom, scales = "free",  space = "free")+
        scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
        theme_bw(tsize) +
        theme(axis.text.x = element_text(face="bold", color="black"),
              axis.text.y = element_text(face="bold", color="black"),
              axis.ticks.y = element_blank(),
              # axis.text.y = element_blank(),
              axis.title.x = element_text(face="bold", color="black"),
              axis.title.y = element_text(face="bold", color="black"),
              strip.text = element_text(face = "bold", color = "black"),
              plot.title = element_text(face="bold"),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        labs(x = "Genomic position (Mb)", y = "") 
    
    # plot genome background
    back <- genodf %>%
        dplyr::filter(sample %in% strainset,
                      chrom == ifelse(chr == "III", "II", "III")) %>%
        dplyr::mutate(bp = end - start) %>%
        dplyr::group_by(sample, gt_name) %>%
        dplyr::summarize(max = sum(bp)) %>%
        dplyr::arrange(desc(max)) %>%
        dplyr::filter(max == first(max)) %>%
        dplyr::mutate(start = 0,
                      end = 15e6,
                      chr = "Genome")  %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sample = factor(sample, levels = strainset)) %>%
        ggplot(.)+
        geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
        facet_grid(~chr, scales = "free",  space = "free")+
        scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
        theme_bw(tsize) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none",
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              strip.text.x = element_text(face = "bold", color = "black"),
              strip.text.y = element_blank())+
        labs(x = "Genomic position (Mb)", y = "")
    
    # return objects
    return(list(chrgeno, back, pheno))
}

# Look for genes in interval
query_genes <- function(region, GO = NULL, strain = "CB4856") {
    
    # filter eqtl to > 5% VE
    eqtlmap2 <- eqtlmap %>%
        dplyr::filter(var_exp >= 0.05)
    
    # how many genes are in the interval?
    all_genes <- cegwas2::query_vcf(region, impact = c("LOW", "MODERATE", "HIGH", "MODIFIER"), samples = strain)
    print(glue::glue("There are {length(unique(all_genes$gene_id))} genes in the interval {region}"))
    
    # how many eQTL map to this region?
    chrom <- stringr::str_split_fixed(region, ":", 2)[,1]
    left_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,1])
    right_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,2])
    
    all_eQTL <- eqtlmap2 %>%
        dplyr::filter(chr == chrom,
                      ci_l_pos < right_pos,
                      ci_r_pos > left_pos)
    print(glue::glue("There are {nrow(all_eQTL)} eQTL ({length(unique(all_eQTL$trait))} traits) that map to {region}"))
    
    # all eQTL probes
    all_eQTL_probes <- eqtl_probes %>%
        dplyr::filter(probe %in% all_eQTL$trait) %>%
        dplyr::left_join(gene_annotations, by = "gene_id")
    
    # which of the eQTL are overlapping with genes in interval?
    eQTL_outside_CI <- all_eQTL_probes %>%
        dplyr::filter(!wbgene %in% all_genes$gene_id)
    print(glue::glue("There are {nrow(all_eQTL)-length(unique(eQTL_outside_CI$wbgene))} genes in the region with an eQTL and {length(unique(eQTL_outside_CI$wbgene))} genes outside the region with an eQTL"))
    
    # Total genes of interest:
    print(glue::glue("There are {length(unique(all_genes$gene_id)) + length(unique(eQTL_outside_CI$wbgene))} total genes of interest."))
    
    # how many of the genes in interval have variation?
    vars <- all_genes %>%
        dplyr::mutate(GT = ifelse(a1 == REF, "ref", "alt")) %>%
        dplyr::filter(GT == "alt")
    
    # genes with protein coding vars
    proteincode <- vars %>%
        dplyr::filter(impact %in% c("MODERATE", "HIGH"))
    print(glue::glue("There are {length(unique(vars$gene_id))}/{length(unique(all_genes$gene_id))} genes in interval with genetic variation, {length(unique(proteincode$gene_id))}/{length(unique(vars$gene_id))} have protein-coding variation"))
    
    # should I look at GO annotations?
    if(!is.null(GO)) {
        # total genes with GO annotations
        go_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_genes$wbgene))}/{length(unique(all_genes$gene_id)) + length(unique(eQTL_outside_CI$wbgene))} genes with {GO} annotation"))
        
        # genes with GO annotations and variation
        go_var <- gene_annotations %>%
            dplyr::filter(wbgene %in% vars$gene_id) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_var$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND genetic variation"))
        
        # genes with GO annotation and protein-coding variation
        go_pcvar <- gene_annotations %>%
            dplyr::filter(wbgene %in% proteincode$gene_id) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_pcvar$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND protein-coding genetic variation"))
        
        # genes with GO annotation and eQTL
        go_eqtl <- gene_annotations %>%
            dplyr::filter(wbgene %in% all_eQTL_probes$wbgene) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_eqtl$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND eQTL"))
        
        # return final dataframe with all info (might be off, only has 133 instead of 134?)
        total_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::mutate(inside_CI = ifelse(wbgene %in% all_genes$gene_id, T, F),
                          eqtl = ifelse(wbgene %in% all_eQTL_probes$wbgene, T, F),
                          vars = ifelse(wbgene %in% vars$gene_id, T, F),
                          pc_vars = ifelse(wbgene %in% proteincode$gene_id, T, F),
                          go_annotation = ifelse(wbgene %in% go_genes$wbgene, T, F))
    } else {
        
        # return final dataframe with all info (might be off, only has 133 instead of 134?)
        total_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::mutate(inside_CI = ifelse(wbgene %in% all_genes$gene_id, T, F),
                          eqtl = ifelse(wbgene %in% all_eQTL_probes$wbgene, T, F),
                          vars = ifelse(wbgene %in% vars$gene_id, T, F),
                          pc_vars = ifelse(wbgene %in% proteincode$gene_id, T, F),
                          go_annotation = NA)
    }
    
    return(total_genes)
}    


##########################################
#           Figure 1                     #
#      RIAILs, linkage, PxG              #
##########################################

riailpheno <- read.csv("data/S3_File.csv")
linkagemapping::load_cross_obj("N2xCB4856cross_full2")
zincmap <- read.csv("data/S4_File.csv")

# trait
t <- "median.EXT"

### plot riail phenos
pheno <- riailpheno %>%
    dplyr::filter(set == 2 | is.na(set), trait == t) %>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain == "N2" ~ "n2",
                                                 strain == "CB4856" ~ "cb",
                                                 TRUE ~ "RIL")) %>%
    dplyr::group_by(strain, trait) %>%
    dplyr::mutate(avg_phen = mean(phenotype, na.rm = T)) %>%
    dplyr::distinct(strain, trait, .keep_all = T) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(norm_pheno = ((avg_phen - min(avg_phen, na.rm = T)) / (max(avg_phen, na.rm = T) - min(avg_phen, na.rm = T)))) %>%
    dplyr::arrange(norm_pheno)
    
pheno$strain <- factor(pheno$strain, levels = unique(pheno$strain))

# number of strains more resistant than N2
n2 <- pheno %>%
    dplyr::filter(strain == "N2") %>%
    dplyr::pull(norm_pheno)
hyperres <- pheno %>%
    dplyr::filter(norm_pheno > n2) %>%
    dplyr::pull(strain) 
length(hyperres)/nrow(pheno) # 13.3%

# number of strains more sensitive than CB
cb <- pheno %>%
    dplyr::filter(strain == "CB4856") %>%
    dplyr::pull(norm_pheno)
hypersen <- pheno %>%
    dplyr::filter(norm_pheno < cb) %>%
    dplyr::pull(strain)
length(hypersen)/nrow(pheno) # 20.3%

(length(hyperres)+length(hypersen))/nrow(pheno) #33.7%

riail <- pheno %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = strain, y = norm_pheno, fill = strain_fill, color = strain_fill) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
    ggplot2::scale_color_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
    theme_bw(tsize) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(face = "bold", color = "black"),
          axis.text.y = element_text(face = "bold", color = "black"),
          legend.position = "none") +
    labs(x = "Strain", y = "Median optical density")

### plot lod

# plot LOD
test <- zincmap %>%
    dplyr::mutate(condtrt = trait,
                  trait = stringr::str_split_fixed(trait, "_", 2)[,2]) %>%
    dplyr::filter(trait == t)
lod <- plot_lods(test)

# add up variance explained
sum(test$var_exp, na.rm = T) # 40.5%

# relative riail pheno
rph <- riailpheno %>% 
    dplyr::filter(trait == t) %>%
    dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = rel_pheno)

# drugcross
drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, rph, set = 2)

# plot pxg
pxg <- plot_pxg(drugcross, test, yaxis = "Median optical density")

# cowplot
cowplot::plot_grid(riail, lod, pxg, nrow = 3, align = "v", axis = "lr", labels = c("A", "B", "C"))

ggsave("figures/Fig1_riail_lod.png", width = 7.5, height = 8)


##########################################
#           Figure 2                     #
#          Chr V NILs                    #
##########################################

zincmap <- read.csv("data/S4_File.csv")
chrVbreakup <- read.csv("data/S11_File.csv")
HTA_stats <- read.csv("data/S10_File.csv")
nil_genotypes <- read.csv("data/S8_File.csv")

# qtl peaks
peaks <- zincmap %>%
    dplyr::filter(trait == "zinc_median.EXT") %>%
    na.omit() %>%
    dplyr::filter(chr == "V")

# plot
trt <- "median.EXT"
tsize <- 12
strainset <- rev(c("N2", "CB4856", "ECA481", "ECA437", "ECA411"))

# regress
regressed <- easysorter::regress(chrVbreakup) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2") ~ "N2",
                                                 strain %in% c("CB4856") ~ "CB",
                                                 TRUE ~ "NIL"),
                  groups = dplyr::case_when(strain == "N2" ~ "N2",
                                            TRUE ~ "CB")) %>%
    dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = rel_pheno)

# calculate NIL effect sizes
test <- regressed %>%
    dplyr::group_by(strain) %>%
    dplyr::summarise(phen = mean(phenotype)) %>%
    tidyr::spread(strain, phen)

(test$ECA411 - test$CB4856) / (test$N2 - test$CB4856)
(test$ECA437 - test$CB4856) / (test$N2 - test$CB4856)
(test$ECA481 - test$CB4856) / (test$N2 - test$CB4856)

am <- lm(phenotype ~ strain, data = regressed) # fix CB4856, N2 as 47.4% effect, ECA411 = 3%, ECA437 = 24%, ECA481 = 32% (diff between ECA481 and ECA411 is 8%, so ECA437 is major locus?)
am

# stats
stats <- HTA_stats %>%
    dplyr::filter(grepl("CB4856", comparison),
                  comparison != "N2-CB4856",
                  experiment == "chrV_breakup",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = stringr::str_split_fixed(comparison, "-", 2)[,1],
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"))

# stats2
stats2 <- HTA_stats %>%
    dplyr::filter(comparison %in% c("ECA481-ECA437", "ECA437-ECA411", "ECA481-ECA411"),
                  experiment == "chrV_breakup",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain1 = stringr::str_split_fixed(comparison, "-", 2)[,1],
                  strain2 = stringr::str_split_fixed(comparison, "-", 2)[,2],
                  yval = dplyr::case_when(comparison == "ECA481-ECA437" ~ 1,
                                          comparison == "ECA437-ECA411" ~ 0.9,
                                          comparison == "ECA481-ECA411" ~ 1.1),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"),
                  movex = dplyr::case_when(comparison == "ECA481-ECA437" ~ -3,
                                           comparison == "ECA437-ECA411" ~ -0.5,
                                           comparison == "ECA481-ECA411" ~ -1.75))

# plot
# plot phenotypes
pheno <- regressed %>%
    dplyr::left_join(stats) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::group_by(strain) %>%
    dplyr::mutate(phen = max(phenotype) + 0.1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(strain = factor(strain, 
                                  levels = strainset),
                  strain_fill = ifelse(strain == "N2", "N2", ifelse(strain == "CB4856", "CB", "NIL")),
                  groups = ifelse(strain == "N2", "N2", "CB")) %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = strain, y = phenotype, fill = strain_fill) +
    ggplot2::geom_jitter(width = 0.1, size = 0.05) +
    ggplot2::geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.2) +
    ggplot2::geom_text(aes(label = sig, y = phen, color = groups), size = tsize/4, angle = -90) +
    ggplot2::geom_segment(data = stats2, aes(x = strain1, xend = strain2, y = yval, yend = yval), inherit.aes = F) +
    ggplot2::geom_text(data = stats2, aes(x = strain1, label = sig, y = yval + 0.02, hjust = movex), angle = -90, size = tsize/4, inherit.aes = F) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
    ggplot2::scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
    ggplot2::theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
           axis.title.x = element_text(face="bold", color="black"),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           strip.text = element_text(face = "bold", color = "black"),
           legend.position = "none",
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank()) +
    coord_flip() +
    labs(x = " ", y = "Median optical density")  +
    facet_grid(~condition)

# plot chrV genotype
chrgeno <- nil_genotypes %>%
    dplyr::filter(chrom == "V",
                  sample %in% strainset) %>%
    dplyr::mutate(chrom = paste0("chr", "V"),
                  sample = factor(sample, levels = strainset)) %>%
    dplyr::distinct() %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(~chrom, scales = "free",  space = "free")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
          axis.text.y = element_text(face="bold", color="black"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(face="bold", color="black"),
          axis.title.y = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black"),
          plot.title = element_text(face="bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(x = "Genomic position (Mb)", y = "") 

# plot genome background
back <- nil_genotypes %>%
    dplyr::filter(sample %in% strainset,
                  chrom == "III") %>%
    dplyr::mutate(bp = end - start) %>%
    dplyr::group_by(sample, gt_name) %>%
    dplyr::summarize(max = sum(bp)) %>%
    dplyr::arrange(desc(max)) %>%
    dplyr::filter(max == first(max)) %>%
    dplyr::mutate(start = 0,
                  end = 15e6,
                  chr = "Genome")  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample = factor(sample, levels = strainset)) %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(~chr, scales = "free",  space = "free")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.text.x = element_text(face = "bold", color = "black"),
          strip.text.y = element_blank())+
    labs(x = "Genomic position (Mb)", y = "")

# plot all together
cowplot::plot_grid(chrgeno +
                       geom_vline(data = peaks, aes(xintercept = pos/1e6), inherit.aes = F) +  # peak marker
                       geom_vline(data = peaks, aes(xintercept = ci_l_pos/1e6), linetype = "dashed") +
                       geom_vline(data = peaks, aes(xintercept = ci_r_pos/1e6), linetype = "dashed"),  # CI
                   back,
                   pheno,
                   nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), align = "h", axis = "bt", labels = c("A", "", "B"))

# save plot
ggsave("figures/Fig2_chrVnil.png", width = 7.5, height = 3)

########## Analysis #############

# analyze genes on chromosome V with zinc annotation
# filtered eqtl with var_exp > 5%
# chrV_genes <- query_genes("V:9620518-11345444", "zinc")


# QTL in ECA437 - MAIN QTL (larger effect size)
chrV2_genes <- query_genes("V:10511995-11345444", "zinc")
test2 <- chrV2_genes %>% dplyr::distinct(wbgene, .keep_all = T)

# QTL in ECA481
chrV1_genes <- query_genes("V:9620518-10511995", "zinc")

##########################################
#           Figure 3                     #
#         dominance HTA                  #
##########################################

dominance_pruned <- read.csv("data/S13_File.csv")
HTA_stats <- read.csv("data/S10_File.csv")
nil_genotypes <- read.csv("data/S8_File.csv")

# plot
trt <- "median.EXT"
tsize <- 12
strainset <- rev(c("C/C", "838/838", "N/838", "859/859","C/859")) 

# regress
regressed <- easysorter::regress(dominance_pruned %>% dplyr::filter(strain %in% strainset)) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = rel_pheno)


# stats
stats <- HTA_stats %>%
    dplyr::filter(comparison %in% c("C/C-859/859", "C/C-C/859"),
                  experiment == "dominance",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = dplyr::case_when(comparison == "N/N-838/838" ~ "838/838",
                                            comparison == "N/N-838/N" ~ "838/N",
                                            comparison == "N/N-N/838" ~ "N/838",
                                            comparison == "C/C-859/859" ~ "859/859",
                                            comparison == "C/C-C/859" ~ "C/859",
                                            TRUE ~ "NA"),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"))

# stats2
stats2 <- HTA_stats %>%
    dplyr::filter(comparison %in% c("N/838-838/838", "C/859-859/859"),
                  experiment == "dominance",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain1 = dplyr::case_when(comparison == "N/838-838/838" ~ "N/838",
                                             comparison == "C/859-859/859" ~ "C/859"),
                  strain2 = dplyr::case_when(comparison == "N/838-838/838" ~ "838/838",
                                             comparison == "C/859-859/859" ~ "859/859"),
                  yval = ifelse(comparison == "N/838-838/838", 1.1, 0.8),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"),
                  movex = ifelse(comparison == "N/838-838/838", 3.5, 2.25))

# plot
# plot phenotypes
pheno <- regressed %>%
    dplyr::left_join(stats) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::group_by(strain) %>%
    dplyr::mutate(phen = max(phenotype) + 0.1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(strain = factor(strain, 
                                  levels = rev(c("N/N", "C/C", "838/838", "838/N", "N/838", "859/859", "859/C", "C/859"))),
                  # strain_fill = ifelse(strain %in% c("N/N", "838/838", "838/N", "N/838"), "N2", "CB")) %>%
                  strain_fill = ifelse(strain == "N/N", "N2", ifelse(strain == "C/C", "CB", "NIL")),
                  groups = ifelse(strain %in% c("N/N", "838/838", "838/N", "N/838"), "N2", "CB")) %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = strain, y = phenotype, fill = strain_fill) +
    ggplot2::geom_jitter(width = 0.1, size = 0.05) +
    ggplot2::geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.2) +
    ggplot2::geom_text(aes(label = sig, y = phen, color = groups), size = tsize/4, angle = -90) +
    ggplot2::geom_segment(data = stats2, aes(x = strain1, xend = strain2, y = yval, yend = yval), inherit.aes = F) +
    ggplot2::geom_text(data = stats2, aes(x = strain1, label = sig, y = yval + 0.1, hjust = movex), size = tsize/4, inherit.aes = F, angle = -90) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
    ggplot2::scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
    ggplot2::theme_bw(tsize) +
    theme(axis.text = element_text(face="bold", color="black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(color = NA),
          strip.text = element_text(face = "bold", color = "black"),
          plot.title = element_text(face="bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          plot.tag = element_text(face = "bold"),
          plot.tag.position = c(-0.0005,1),
          axis.title.x = element_text(face = "bold", color = "black"),
          panel.grid.major = element_blank()) +
    ggplot2::labs(x = "", y = "Median optical density") +
    ggplot2::coord_flip() +
    facet_grid(~condition)


# plot chromosome genotypes
df <- data.frame(strain = strainset) %>%
    dplyr::mutate(sample = strain,
                  sample = gsub("C", "CB4856", sample),
                  sample = gsub("N", "N2", sample),
                  sample = gsub("838", "ECA838", sample),
                  sample = gsub("859", "ECA859", sample)) %>%
    tidyr::separate_rows(sample, sep = "/") %>%
    dplyr::left_join(nil_genotypes) %>%
    dplyr::filter(chrom == "III") %>%
    dplyr::mutate(chrom = "chrIII") %>%
    dplyr::mutate(strain2 = dplyr::case_when(strain == "N/N"~ "N2/\nN2",
                                             strain == "C/C" ~"CB4856/\nCB4856",
                                             strain == "838/838" ~ "ECA838/\nECA838",
                                             strain == "859/859" ~ "ECA859/\nECA859",
                                             strain == "N/838" ~ "N2/\nECA838",
                                             strain == "C/859" ~ "CB4856/\nECA859"))
df <- cbind(df, ystrain = c(2, 1, 1, 3, 3, 4, 4, 6, 5, 5, 5, 5, 7, 7, 7, 7, 8, 8, 8, 8, 9, 10))


chrIII_plot <- df %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = as.factor(ystrain), xend = end/1e6, yend = as.factor(ystrain), color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(factor(strain2, levels = c("N2/\nN2", "CB4856/\nCB4856", "ECA838/\nECA838", "N2/\nECA838", "ECA859/\nECA859", "CB4856/\nECA859"))~chrom, scales = "free",  space = "free", switch = "y")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(face="bold", color="black"),
          axis.title.y = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black"),
          plot.title = element_text(face="bold"),
          legend.position = "none",
          plot.tag = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.background.y = element_rect(fill = "white", color = "NA"),
          strip.text.y = element_text(angle = 180),
          panel.grid.major = element_blank())+
    labs(x = "Genomic position (Mb)", y = "") 

# DREADED GROBS
elements <- ggplot2::ggplotGrob(chrIII_plot)
panels = subset(elements$layout, grepl("axis", elements$layout$name), t:r)
chrIII_plot2 <- gtable::gtable_add_grob(elements, 
                                        grid::rectGrob(gp = grid::gpar(col = "black", lwd = 1, fill = NA)),
                                        t = 8, l = 6, b = 16, r = 7)
chrIII_plot3 <- ggplotify::as.ggplot(chrIII_plot2)

# background genome
genome_plot <- data.frame(strain = strainset) %>%
    dplyr::mutate(sample = dplyr::case_when(strain %in% c("C/859", "859/859", "C/C") ~ "CB4856,CB4856",
                                            TRUE ~ "N2,N2")) %>%
    tidyr::separate_rows(sample, sep = ",") %>%
    dplyr::left_join(nil_genotypes) %>%
    dplyr::filter(chrom == "I") %>%
    dplyr::mutate(chrom = "Genome") %>%
    tibble::rownames_to_column(var = "ystrain") %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = as.factor(ystrain), xend = end/1e6, yend = as.factor(ystrain), color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(factor(strain, levels = rev(strainset))~chrom, scales = "free",  space = "free", switch = "y")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text.x = element_text(color = NA),
          axis.title.x = element_text(color = NA),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_line(color = NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(face = "bold", color = "black"),
          strip.text.y = element_blank())+
    labs(x = "Genomic position (Mb)", y = "")

# GROBS
elements <- ggplot2::ggplotGrob(genome_plot)
panels = subset(elements$layout, grepl("axis", elements$layout$name), t:r)
genome_plot2 <- gtable::gtable_add_grob(elements, 
                                        grid::rectGrob(gp = grid::gpar(col = "black", lwd = 1, fill = NA)),
                                        t = 8, l = 6, b = 16, r = 7)
genome_plot3 <- ggplotify::as.ggplot(genome_plot2)


test1 <- cowplot::plot_grid(chrIII_plot3, genome_plot3, ncol = 2, rel_widths = c(3.25, 1), labels = c("A", ""))

cowplot::plot_grid(test1, pheno, ncol = 2, rel_widths = c(4, 4.5), labels = c("", "B"))

ggsave("figures/Fig3_dominance.png", width = 7.5, height = 5)



########## Analysis #############

# analyze genes on chromosome III with zinc annotation
chrIII_genes <- query_genes("III:4664-597553", "zinc")


##########################################
#           Figure 4                     #
#          sqst5 eQTL                    #
##########################################

# ver-2 probe: A_12_P104472
sqst5_mediation <- read.csv("data/S16_File.csv")
zincmap <- read.csv("data/S4_File.csv")
riailpheno <- read.csv("data/S3_File.csv")
sqst5_map <- read.csv("data/S15_File.csv")
linkagemapping::load_cross_obj("N2xCB4856cross_full2")

# get set2 mapping
map <- zincmap %>%
    na.omit() %>%
    dplyr::filter(trait == "zinc_median.EXT", 
                  chr == "III") 

med <- sqst5_mediation %>%
    dplyr::mutate(abs_est = abs(estimate)) %>%
    tidyr::separate(marker, c("chrom", "pos"), "_") %>%
    dplyr::mutate(pos = as.numeric(pos)) %>%
    dplyr::filter(var == "med") %>%
    dplyr::arrange(pos) %>%
    dplyr::mutate(sqst5 = ifelse(probe == "A_12_P104472", "yes", "no"),
                  condition = stringr::str_split_fixed(trait, "_", 2)[,1]) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(q90 = quantile(abs_est, probs = 0.9)[[1]],
                  q95 = quantile(abs_est, probs = 0.95)[[1]],
                  q99 = quantile(abs_est, probs = 0.99)[[1]]) %>%
    ggplot(.) +
    aes(x = pos/1e6, y = abs_est, fill = sqst5, shape = sqst5, size = sqst5) +
    geom_point(aes(alpha = prob), color = "black") +
    scale_alpha_continuous(range = c(1, 0.1)) +
    theme_bw(tsize) +
    scale_fill_manual(values = c("yes" = "red", "no" = "black")) +
    scale_shape_manual(values = c("yes" = 23, "no" = 21)) +
    scale_size_manual(values = c("yes" = 3, "no" = 1)) +
    labs(x = "Genomic position (Mb)", y = "Mediation estimate") +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.text = element_text(face="bold", color="black"),
          axis.title = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black")) +
    geom_hline(aes(yintercept = q95), color = "grey") +
    facet_wrap(~condition, scales = "free", nrow = 2)

# sqst-5 eQTL LOD
# load plot_lods() and plot_pxg() functions from above
lod <- plot_lods(sqst5_map %>%
              dplyr::mutate(condtrt = trait,
                            trait = "sqst-5"))

# cross obj
drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, eqtlpheno %>%
                                            dplyr::filter(probe == "A_12_P104472") %>%
                                            dplyr::mutate(trait = probe,
                                                          condition = "expression",
                                                          phenotype = expression), set = 1)
# plot pxg
pxg <- plot_pxg(drugcross, sqst5_map %>%
             dplyr::mutate(trait = "expression_A_12_P104472",
                           condtrt = trait,
                           condition = "expression"),
             yaxis = "Expression")


bcd <- cowplot::plot_grid(pxg, med, nrow = 1, align = "h", axis = "b", labels = c("B", "C"))

cowplot::plot_grid(lod, bcd, nrow = 2, align = "v", labels = c("A", ""))

ggsave("figures/Fig4_eqtlmediation.png", width = 7.5, height = 5)


##########################################
#           Figure 5                     #
#        sqst5 hemi HTA                  #
##########################################

# load file
sqst5_hemi_pruned <- read.csv("data/S19_File.csv")
HTA_stats <- read.csv("data/S10_File.csv")
nil_genotypes <- read.csv("data/S8_File.csv")

# trait
trt <- "median.EXT"

strainset <- rev(unique(regressed$strain))

# regress
regressed <- easysorter::regress(sqst5_hemi_pruned %>%
                                     dplyr::filter(strain %in% c("N2/N2", "CB4856/CB4856", "ECA859/ECA859", "ECA2517/ECA2517", "ECA859/ECA2517"))) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = rel_pheno)

# stats comparison
all_comparisons <- HTA_stats %>%
    dplyr::filter(experiment == "sqst5_hemi") %>%
    dplyr::pull(comparison) %>%
    unique()

# comparisons to keep
cb_compare <- grep("CB4856/CB4856", all_comparisons, value = T)
eca859_compare <- grep("ECA859/ECA859", all_comparisons, value = T)

# stats
stats <- HTA_stats %>%
    dplyr::filter(comparison %in% c(cb_compare, eca859_compare),
                  !grepl("N2/N2", comparison),
                  experiment == "sqst5_hemi",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = dplyr::case_when(grepl("CB4856/CB4856", comparison) ~ stringr::str_split_fixed(comparison, "-", 2)[,1],
                                            #grepl("ECA859/ECA859", comparison) ~ stringr::str_split_fixed(comparison, "-", 2)[,2],
                                            TRUE ~ "NA"),
                  group = dplyr::case_when(grepl("CB4856/CB4856", comparison) ~ "cb",
                                           grepl("ECA859/ECA859", comparison) ~ "859",
                                           TRUE ~ "NA"),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"))

# stats2
stats2 <- HTA_stats %>%
    dplyr::filter(comparison %in% c("ECA859/ECA859-ECA2517/ECA2517", "ECA859/ECA859-ECA859/ECA2517", "ECA859/ECA2517-ECA2517/ECA2517"),
                  experiment == "sqst5_hemi",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain1 = stringr::str_split_fixed(comparison, "-", 2)[,1],
                  strain2 = stringr::str_split_fixed(comparison, "-", 2)[,2],
                  yval = dplyr::case_when(comparison == "ECA859/ECA859-ECA2517/ECA2517" ~ 0.9,
                                          comparison == "ECA859/ECA2517-ECA2517/ECA2517" ~ 1.1,
                                          comparison == "ECA859/ECA859-ECA859/ECA2517" ~ 1.3),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"),
                  movex = dplyr::case_when(comparison == "ECA859/ECA859-ECA2517/ECA2517" ~ -2.25,
                                           comparison == "ECA859/ECA2517-ECA2517/ECA2517" ~ 4,
                                           comparison == "ECA859/ECA859-ECA859/ECA2517" ~ -4))

# plot phenotypes
pheno <- regressed %>%
    dplyr::filter(trait == trt, condition == "zinc") %>%
    dplyr::group_by(strain, condition) %>%
    dplyr::mutate(phen = max(phenotype) + 0.1) %>%
    dplyr::ungroup() %>%
    dplyr::full_join(stats, by = c("strain", "condition", "trait")) %>%
    dplyr::mutate(strain = factor(strain, 
                                  # levels = c("N2/N2", "CB4856/CB4856", "ECA859/ECA859", "ECA2517/ECA2517", "ECA859/CB4856", "ECA2517/CB4856", "ECA859/ECA2517")),
                                  levels = rev(c("N2/N2", "CB4856/CB4856", "ECA859/ECA859", "ECA2517/ECA2517", "ECA859/ECA2517"))),
                  # levels = c("N2/N2", "CB4856/CB4856", "ECA859/ECA859", "ECA2517/ECA2517", "ECA859/CB4856", "ECA2517/CB4856", "CB4856/ECA2517", "ECA859/ECA2517", "ECA2517/ECA859")),
                  strain_fill = ifelse(strain == "N2/N2", "n2", ifelse(strain == "CB4856/CB4856", "cb", "nil"))) %>%
    dplyr::filter(!is.na(strain)) %>%
    # dplyr::mutate(strain = gsub("/", "/\n", strain)) %>%
    ggplot(.) +
    aes(x = strain, y = phenotype, fill = strain_fill) +
    geom_jitter(width = 0.1, size = 0.05) +
    geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.2) +
    ggplot2::geom_text(aes(label = sig, y = phen, color = group), size = tsize/4, angle = -90) +
    ggplot2::geom_segment(data = stats2, aes(x = strain1, xend = strain2, y = yval, yend = yval), inherit.aes = F) +
    ggplot2::geom_text(data = stats2, aes(x = strain1, label = sig, y = yval + 0.1, hjust = movex), size = tsize/4, inherit.aes = F, angle = -90) +
    scale_fill_manual(values = c("n2" = "orange", "cb" = "blue", "nil" = "grey")) +
    scale_color_manual(values = c("859" = "orange", "cb" = "blue")) +
    theme_bw(tsize) +
    scale_x_discrete(labels = c("N2/\nN2", "CB4856/\nCB4856", "ECA859/\nECA859",
                                "ECA2517/\nECA2517", "ECA859/\nECA2517")) +
    theme(axis.text = element_text(face="bold", color="black"),
          axis.title.x = element_text(face="bold", color="black"),
          axis.title.y = element_text(color = NA),
          axis.text.y = element_blank(),
          strip.text = element_text(face = "bold", color = "black"),
          legend.position = "none",
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    labs(x = " ", y = "Median optical density")+
    ggplot2::coord_flip() +
    facet_grid(~condition)

# pheno

# add fake sqst del genotypes
eca2517 <- nil_genotypes %>%
    dplyr::filter(sample == "ECA859") %>%
    dplyr::mutate(sample = "ECA2517")
eca2518 <- nil_genotypes %>%
    dplyr::filter(sample == "ECA859") %>%
    dplyr::mutate(sample = "ECA2518")
nil_genotypes <- nil_genotypes %>%
    dplyr::bind_rows(eca2517, eca2518)

# plot chromosome genotypes
df <- data.frame(strain = strainset) %>%
    dplyr::mutate(sample = strain) %>%
    tidyr::separate_rows(sample, sep = "/") %>%
    dplyr::left_join(nil_genotypes) %>%
    dplyr::filter(chrom == "III") %>%
    dplyr::mutate(chrom = "chrIII") %>%
    dplyr::mutate(strain2 = dplyr::case_when(strain == "N2/N2"~ "N2/\nN2",
                                             strain == "CB4856/CB4856" ~"CB4856/\nCB4856",
                                             strain == "ECA859/ECA859" ~ "ECA859/\nECA859",
                                             strain == "ECA2517/ECA2517" ~ "ECA2517/\nECA2517",
                                             strain == "ECA859/ECA2517" ~ "ECA859/\nECA2517"))
df <- cbind(df, ystrain = c(2,2,1,1,3,3,4,4,5,5,6,6,7,8,9,10))

chrIII_plot <- df %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = as.factor(ystrain), xend = end/1e6, yend = as.factor(ystrain), color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(factor(strain2, levels = c("N2/\nN2", "CB4856/\nCB4856", "ECA859/\nECA859", "ECA2517/\nECA2517", "ECA859/\nECA2517"))~chrom, scales = "free",  space = "free", switch = "y")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(face="bold", color="black"),
          axis.title.y = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black"),
          plot.title = element_text(face="bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.background.y = element_rect(fill = "white", color = "NA"),
          strip.text.y = element_text(angle = 180),
          panel.grid.major = element_blank())+
    labs(x = "Genomic position (Mb)", y = "") +
    geom_point(data = data.frame(sample = c(rep("ECA2517", 3)), ystrain = factor(c(1, 3, 4)), strain2 = c("ECA859/\nECA2517", "ECA2517/\nECA2517", "ECA2517/\nECA2517")),
               aes(x = 0.146, y = ystrain), shape = 24, color = "black", fill = "grey", size = 2)

# DREADED GROBS
elements <- ggplot2::ggplotGrob(chrIII_plot)
panels = subset(elements$layout, grepl("axis", elements$layout$name), t:r)
chrIII_plot2 <- gtable::gtable_add_grob(elements, 
                                        grid::rectGrob(gp = grid::gpar(col = "black", lwd = 1, fill = NA)),
                                        t = 8, l = 6, b = 16, r = 7)
chrIII_plot3 <- ggplotify::as.ggplot(chrIII_plot2)


# background genome
genome_plot <- data.frame(strain = strainset) %>%
    dplyr::mutate(sample = dplyr::case_when(strain == "N2/N2" ~ "N2,N2", 
                                            TRUE ~ "CB4856,CB4856")) %>%
    tidyr::separate_rows(sample, sep = ",") %>%
    dplyr::left_join(nil_genotypes) %>%
    dplyr::filter(chrom == "I") %>%
    dplyr::mutate(chrom = "Genome") %>%
    tibble::rownames_to_column(var = "ystrain") %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = as.factor(ystrain), xend = end/1e6, yend = as.factor(ystrain), color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(factor(strain, levels = rev(strainset))~chrom, scales = "free",  space = "free", switch = "y")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text.x = element_text(color = NA),
          axis.title.x = element_text(color = NA),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_line(color = NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(face = "bold", color = "black"),
          strip.text.y = element_blank())+
    labs(x = "Genomic position (Mb)", y = "")

# GROBS
elements <- ggplot2::ggplotGrob(genome_plot)
panels = subset(elements$layout, grepl("axis", elements$layout$name), t:r)
genome_plot2 <- gtable::gtable_add_grob(elements, 
                                        grid::rectGrob(gp = grid::gpar(col = "black", lwd = 1, fill = NA)),
                                        t = 8, l = 6, b = 16, r = 7)
genome_plot3 <- ggplotify::as.ggplot(genome_plot2)

test1 <- cowplot::plot_grid(chrIII_plot3, genome_plot3, ncol = 2, rel_widths = c(3.3, 1), labels = c("A", ""))

cowplot::plot_grid(test1, pheno, ncol = 2, rel_widths = c(4, 4), labels = c("", "B"))
ggsave("figures/Fig5_sqsthemi.png", width = 7.5, height = 5)


##########################################
#           Figure 6                     #
#          SQST5 variant tree            #
##########################################

WI_sqst5_sv <- read.csv("data/S23_File.csv")
tree <- ape::read.tree("~/Dropbox/AndersenLab/LabFolders/Katie/projects/zinc/manuscript/data/processed/FileS24_nj_sqst.tree")

# gene model
sqst5_N2 <- data.frame(starts = c(147608,147314), ends = c(147361, 146983), strain = "N2", txstart = 147608, txend = 146983)

introns_N2 <- data.frame(starts = c(147361), ends = c(147314)) %>%
    dplyr::mutate(midpoint = ((ends - starts) / 2) + starts)

endUTR_N2 <- data.frame(x = c(146917, 146982, 146982), y = c(0, 1, -1))
begUTR_N2 <- data.frame(starts = c(147620), ends = c(147608), strain = "N2")

gene_model <- ggplot2::ggplot(sqst5_N2)+
    # coding region
    # ggplot2::geom_segment(aes(x = txend, xend = txstart, y = 0, yend = 0), color = "black")+
    # exons
    ggplot2::geom_rect(aes(xmin =  starts, xmax = ends, ymin = -1 , ymax = 1), color = "black", alpha = 0.5, fill = "grey70")+
    # introns
    ggplot2::geom_segment(aes(x = starts, y = 1, xend = midpoint, yend = 2), data = introns_N2, color = "black")+
    ggplot2::geom_segment(aes(x = midpoint, y = 2, xend = ends, yend = 1), data = introns_N2, color = "black")+
    # upstream UTR
    ggplot2::geom_rect(aes(xmin = starts, xmax = ends, ymin = -1, ymax = 1), data = begUTR_N2, fill= "white", color = "black")+
    # downstream UTR
    ggplot2::geom_polygon(aes(x = x, y = y), data = endUTR_N2, fill = "gray60", color = "black")+
    ggplot2::theme_void() +
    # zz type domain
    ggplot2::geom_rect(aes(xmin = 147204, xmax = 147064, ymin = -1, ymax = 1), fill = "darkgreen", color = "black", alpha = 0.5) +
    ggplot2::geom_text(label = "ZZ", color = "white", aes(x = 147134, y = 0), size = 5) +
    # deletion in CB4856
    ggplot2::geom_rect(aes(xmin = 147186, xmax = 147076, ymin = -1.25, ymax = -2), fill = "navy", color = "black", size = 0.3) +
    ggplot2::geom_text(label = "CB4856 del", color = "white", fontface = "bold", aes(x = 147131, y = -1.625), size = 3) 
# CYS-2X-CYS residues
# ggplot2::geom_text(label = "C", color = "white", aes(x = c(147185,147176,147140,147131,147113,147104), y = c(1,1,1,1,1,1)), size = 5)

# tree
# strains with deletion:
delstrains <- WI_sqst5_sv %>%
    dplyr::filter(CB4856_del == T) %>%
    dplyr::pull(strain)

otherdel <- WI_sqst5_sv %>%
    dplyr::filter(CB4856_del == F, sqst5_sv == T) %>%
    dplyr::pull(strain)

missing <- WI_sqst5_sv %>%
    dplyr::filter(is.na(sqst5_sv))%>%
    dplyr::pull(strain)

treeplot <- as.tibble(phytools::midpoint.root(tree)) %>% 
    dplyr::rename(strain=label) %>%
    dplyr::mutate(label = dplyr::case_when(strain == "N2" ~ "N2",
                                           strain == "CB4856" ~ "CB",
                                           strain %in% delstrains ~ "del",
                                           strain %in% otherdel ~ "other",
                                           strain %in% missing ~ "missing",
                                           TRUE ~ "nodel")) %>%
    dplyr::select(parent,node,branch.length,label, strain) %>%
    tidytree::as.phylo(tbl_tree_admix_KSE, use.labels = T) %>%
    ggtree::ggtree(size = 0.3,color="gray51") +
    ggtree::geom_tippoint(aes(fill=label), 
                          shape=21, 
                          size=2)  +  
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "navy", "del" = "mediumvioletred", 
                                          "nodel" = "grey", "other" = "pink", "missing" = "white"),
                               labels = c("N2" = "N2", "CB" = "CB4856", "del" = "Deletion", 
                                          "nodel" = "Wild-type", "other" = "Variation", "missing" = "Missing"),
                               name = expression(paste(italic("sqst-5"), " deletion"))) +
    ggplot2::scale_x_reverse() +
    ggplot2::theme(legend.position = c(0.8,0.7))
# treeplot

cowplot::plot_grid(cowplot::plot_grid(gene_model, NULL, nrow = 1, rel_widths = c(3, 1)), treeplot, nrow = 2, rel_heights = c(1, 4), labels = c("A", "B"))
ggsave("figures/Fig6_phylo.png", height = 5, width = 7.5)


##########################################
#           Figure 7                     #
#             GWAS                       #
##########################################

# make plots like linkage
WI_sqst5_sv <- read.csv("data/S23_File.csv")
gwaspheno <- read.csv("data/S26_File.csv")
gwasmap <- read.csv("data/S27_File.csv")

# plot gwas for median.norm.ext
t <- "median.norm.EXT"

### plot WI phenos
pheno <- gwaspheno %>%
    dplyr::filter(trait == t) %>%
    dplyr::left_join(WI_sqst5_sv, by = "strain") %>%
    tidyr::drop_na(sqst5_sv) %>%
    dplyr::mutate(norm_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::arrange(norm_pheno)

pheno$strain <- factor(pheno$strain, levels = unique(pheno$strain))

wi_pheno <- pheno %>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain == "N2" ~ "N2",
                                                 strain == "CB4856" ~ "CB4856",
                                                 CB4856_del == T ~ "Deletion",
                                                 sqst5_sv == T ~ "Variation",
                                                 TRUE ~ "Wild-type")) %>%
    dplyr::mutate(strain_fill = factor(strain_fill, levels = c("Wild-type", "N2", "CB4856", "Deletion", "Variation"))) %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = strain, y = norm_pheno, fill = strain_fill) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB4856" = "navy", "Wild-type" = "grey", 
                                          "Deletion" = "mediumvioletred", "Variation" = "pink"),
                               guide = guide_legend(direction = 'horizontal')) +
    theme_bw(tsize) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(face = "bold", color = "black"),
          axis.text.y = element_text(face = "bold", color = "black"),
          legend.position = c(0.35, 0.9),
          legend.title = element_blank()) +
    labs(x = "Strain", y = t)

# GWAS mapping
gwas_map <- gwasmap %>%
    dplyr::filter(trait == glue::glue("zinc_{t}") & CHROM != "MtDNA") %>%
    dplyr::distinct(marker, .keep_all = T) %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = POS/1e6, y = log10p) +
    ggplot2::scale_color_manual(values = c("0" = "black", "1" = "red", "2" = "hotpink3")) +
    ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6,    # this is the plot boundary for LD and gene plots
                                    xmax = endPOS/1e6,    # this is the plot boundary for LD and gene plots
                                    ymin = 0, 
                                    ymax = Inf, 
                                    fill = "skyblue"), 
                       color = "blue",fill = "skyblue",linetype = 2, alpha=.3)+
    ggplot2::geom_hline(ggplot2::aes(yintercept = BF), color = "gray", alpha = .75, size = 1) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF), color = "gray", alpha = .75, size = 1, linetype = 2) +
    ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG)) ) +
    ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
    ggplot2::theme_bw(12) +
    ggplot2::theme(
        axis.text = ggplot2::element_text(color = "black", face = "bold"),
        axis.title = ggplot2::element_text(face = "bold", color = "black"),
        strip.text = ggplot2::element_text(face = "bold", color = "black"),
        plot.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
    ggplot2::labs(x = "Genomic position (Mb)",
                  y = expression(-log[10](italic(p))))

# pxg from gwas
peaks <- gwasmap %>%
    dplyr::filter(trait == glue::glue("zinc_{t}") & CHROM != "MtDNA") %>%
    na.omit()

pxg <- peaks %>% # probably the na.omit is not needed any more
    dplyr::mutate(value=as.numeric(value),
                  allele = factor(allele, levels = c(-1,1), labels = c("REF","ALT"))) %>% 
    dplyr::mutate(norm_pheno = ((value - min(value, na.rm = T)) / (max(value, na.rm = T) - min(value, na.rm = T))),
                  value = norm_pheno) %>%
    dplyr::group_by(allele) %>%
    ggplot(.)+
    aes(x = allele, y = value, fill = allele)+
    geom_jitter(size = 0.5, width = 0.1) +
    geom_boxplot(alpha = 0.8, outlier.color = NA) +
    scale_fill_manual(values = c("REF" = "grey", "ALT" = "navy")) +
    ggplot2::facet_grid(~factor(marker, unique(peaks$marker)), scales = "free") + 
    ggplot2::theme_bw(tsize) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", color = "black"), 
                   axis.text.y = ggplot2::element_text(face = "bold", color = "black"), 
                   axis.title.x = ggplot2::element_text(face = "bold", color = "black", vjust = -0.3),
                   axis.title.y = ggplot2::element_text(face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
                   strip.text.y = ggplot2::element_text(face = "bold", color = "black"),
                   plot.title = ggplot2::element_blank(), 
                   legend.position = "none", 
                   panel.grid = element_blank(),
                   panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
    ggplot2::labs(x = "Genotype at QTL", y = t)
    
# cowplot
cowplot::plot_grid(wi_pheno, gwas_map, pxg, nrow = 3, align = "v", axis = "lr", labels = c("A", "B", "C"))
ggsave("figures/Fig7_gwas.png", width = 7.5, height = 8)

