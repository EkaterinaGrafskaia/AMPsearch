################# Physico-chemical properties #################
library("seqinr")
library("protr")
library("Peptides")
library("segmented")
xf <- read.fasta(file = "Proteins.fasta", seqtype = "AA", as.string = TRUE)
pcpFunc <- function(hn) {
  l <- length(hn)
  result <- list()
  for (i in 1:l){
    seqname <- attr(hn[[i]],"name")
    lseq <- lengthpep(hn[[i]]) # integer
    molw <- round(mw(hn[[i]]), digits=2) 
    ain <- round(aIndex(hn[[i]]), digits=2) #numeric
    boman <- round(boman(hn[[i]]), digits=2) #numeric
    gravy <- round(hydrophobicity(hn[[i]],"KyteDoolittle"), digits=2)
    gravya <- round(hydrophobicity(hn[[i]],"Argos"), digits=2)
    gravyh <- round(hydrophobicity(hn[[i]],"HoppWoods"), digits=2)
    gravyj <- round(hydrophobicity(hn[[i]],"Janin"), digits=2)
    helic <- round(hmoment(hn[[i]], angle = 100, window = 11),digits=2)
    beta <- round(hmoment(hn[[i]], angle = 160, window = 11),digits=2)
    piemb <- round(pI(hn[[i]],"EMBOSS"), digits=2) #numeric
    chseq <- round(charge(hn[[i]],pH=seq(from = 5,to = 9,by = 2), pKscale="EMBOSS"), digits=2) # numeric
    aacomp <-  aaComp(hn[[i]])
    dalm <- round(aacomp[[1]][3]*100/lseq, digits=2)
    darm <- round(aacomp[[1]][4]*100/lseq, digits=2)
    dnpm <- round(aacomp[[1]][5]*100/lseq, digits=2)
    dpom <- round(aacomp[[1]][6]*100/lseq, digits=2)
    dchm <- round(aacomp[[1]][7]*100/lseq, digits=2)
    dbm <- round(aacomp[[1]][8]*100/lseq, digits=2)
    dacm <- round(aacomp[[1]][9]*100/lseq, digits=2)
    result[[i]] <- list(seqname, hn[i],lseq, molw, ain, boman,gravy,gravya, gravyh,gravyj,
                        helic,beta,piemb, chseq,  dalm, darm, dnpm, dpom,dchm, dbm, dacm)
  }
  output <- matrix(unlist(result), byrow=TRUE, nrow=length(result) )
  colnames(output) <- c("Seq_name", "Sequence", "Length", "Mol_weight", "Aliphat_index","Boman_index","GRAVY_hydrophobicity_index_KD",
                        "GRAVY_hydrophobicity_index_A","GRAVY_hydrophobicity_index_HW","GRAVY_hydrophobicity_index_J",
                        "Hydrophobic_moment_helical","Hydrophobic_moment_bstrand",
                        "pI_EMBOSS","Charge_pH=5","Charge_pH=7","Charge_pH=9",
                        "Aliphatic_Mole", "Aromatic_Mole", "Non_polar_Mole","Polar_Mole","Charged_Mole","Basic_Mole", "Acidic_Mole")
  write.table(file="Result_peptides_phys_chem_prop.csv", output, sep="\t")
}

pcpFunc(xf)


xf <- readFASTA('Protein.fasta')
aacFunc <- function(hn) {
  l <- length(hn)
  result <- list()
  for (i in 1:l){
    aacomp <-  extractAAC(hn[[i]])
    darg <- 100*round(aacomp[[2]], digits=2)
    dhis <- 100*round(aacomp[[9]], digits=2)
    dleu <- 100*round(aacomp[[11]], digits=2)
    dlys <- 100*round(aacomp[[12]], digits=2)
    dpro <- 100*round(aacomp[[15]], digits=2)
    dtrp <- 100*round(aacomp[[18]], digits=2)
    result[[i]] <- list(hn[i],darg,dhis,dleu,dlys,dpro,dtrp)
  }
  output <- matrix(unlist(result), byrow=TRUE, nrow=length(result) )
  colnames(output) <- c("Sequence","Arg",	"His","Leu","Lys", "Pro",	"Trp")
  write.table(file="Result_aacomp.csv", output, sep="\t")
}
aacFunc(xf)



