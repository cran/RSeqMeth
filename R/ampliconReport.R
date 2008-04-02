ampliconReport <- function (fnames = NA, minMass=1500, maxMass=7000)
{
ampliconReportT <- function (fnames = NA, minMass, maxMass)
{
    generateReport <- function(fname) {
        AMPa <- parseSequence(fname)
        AMPc <- findCpGs(AMPa)
        AMPu <- findUniqueUnits(AMPc, minMass, maxMass, skipOL = TRUE)
        AMPuOL <- findUniqueUnits(AMPc, minMass, maxMass, skipOL = FALSE)
        allUnits <- list(Frags = AMPa, Peaks = makePeaks(AMPa))
        CpGUnits <- list(Frags = AMPc, Peaks = makePeaks(AMPc))
        UniqueUnits <- list(Frags = AMPu, Peaks = makePeaks(AMPc))
        UniqueUnitsOL <- list(Frags = AMPuOL, Peaks = makePeaks(AMPc))
        header <- rep("", 11)
        header[1] = paste("T Cleavage report of amplicon", fname)
        header[2] = paste("No of Fragments:", length(allUnits$Frags))
        header[3] = paste("No of Peaks:", length(allUnits$Peaks))
        CpGCount <- 0
        nCpGUnits <- length(CpGUnits$Frags)
        for (i in 1:length(CpGUnits$Frags)) CpGCount = CpGCount +
            CpGUnits$Frags[[i]]$CpGs
        header[4] = paste(CpGCount, "CpGs", "in", nCpGUnits,
            "CpG Units")
        CpGCount <- 0
        nUniqueUnits <- length(UniqueUnitsOL$Frags)
        for (i in 1:length(UniqueUnitsOL$Frags)) CpGCount = CpGCount +
            UniqueUnitsOL$Frags[[i]]$CpGs
        header[5] = paste(CpGCount, "CpGs", "in", nUniqueUnits,
            "informative CpG Units (including overlaps)")
        header[7] = "Predicted T Cleavage fragments - top (template) and bottom (transcribed) strand"
        header[8:9] = showCleavage(AMPa)
        for (i in nchar(fname):1) if (substr(fname, i, i) ==
            ".")
            break
        fname2 <- paste(substr(fname, 1, (i - 1)), " T Report.csv",
            sep = "")
        options(warn = -1)
        write(header, fname2, append = FALSE)
        nfrags <- length(CpGUnits$Frags)
        cnames <- c("Name", "Position", "Sequence", "CpGs", "Base mass",
            "Desc", "<15.9Da")
        tname <- vector(mode = "character", length = nfrags)
        tPos <- vector(mode = "numeric", length = nfrags)
        tSeq <- vector(mode = "character", length = nfrags)
        tNo <- vector(mode = "numeric", length = nfrags)
        tBase <- vector(mode = "numeric", length = nfrags)
        tDesc <- vector(mode = "character", length = nfrags)
        tClose <- rep(NA, nfrags)
        for (i in 1:nfrags) {
            tname[i] = CpGUnits$Frags[[i]]$Name
            tNo[i] = CpGUnits$Frags[[i]]$CpGs
            tPos[i] = CpGUnits$Frags[[i]]$Pos
            tSeq[i] = CpGUnits$Frags[[i]]$Seq
            tBase[i] = CpGUnits$Frags[[i]]$Mass[1]
            tDesc[i] = CpGUnits$Frags[[i]]$CpGSites
            noPeaks <- length(CpGUnits$Frags[[i]]$Mass)
            Cfrags <- vector(length = noPeaks)
            Afrags <- vector(length = noPeaks)
            for (j in 1:noPeaks) {
                Cfrags[j] = length(findFrag(CpGUnits$Frags[[i]]$Mass[j],
                  CpGUnits$Peaks)$Frag)
                Afrags[j] = length(findFrag(CpGUnits$Frags[[i]]$Mass[j],
                  allUnits$Peaks)$Frag)
                temp <- findFragMultiple((CpGUnits$Frags[[i]]$Mass[j] -
                  8), allUnits$Peaks, tol = 7.999)
                temp <- c(temp, findFragMultiple((CpGUnits$Frags[[i]]$Mass[j] +
                  8), allUnits$Peaks, tol = 7.999))
                if (length(temp) > 0)
                  for (k in 1:length(temp)) {
                    dist <- abs(CpGUnits$Frags[[i]]$Mass[j] -
                      temp[[k]]$Mass)
                    if (is.na(tClose[i]) || (tClose[i] > dist))
                      tClose[i] = dist
                  }
            }
        }
        fTable <- data.frame(tname, tPos, tSeq, tNo, tBase, tDesc,
            tClose)
        names(fTable) <- cnames
        write("All CpG_Units", fname2, append = TRUE)
        write.csv(fTable, fname2, append = TRUE)
        nfrags <- length(UniqueUnitsOL$Frags)
        cnames <- c("Name", "Base mass", "Desc", "<15.9Da")
        tname <- vector(mode = "character", length = nfrags)
        tBase <- vector(mode = "numeric", length = nfrags)
        tDesc <- vector(mode = "character", length = nfrags)
        tClose <- rep(NA, nfrags)
        for (i in 1:nfrags) {
            tname[i] = UniqueUnitsOL$Frags[[i]]$Name
            tBase[i] = UniqueUnitsOL$Frags[[i]]$Mass[1]
            tDesc[i] = UniqueUnitsOL$Frags[[i]]$CpGSites
            noPeaks <- length(UniqueUnitsOL$Frags[[i]]$Mass)
            Cfrags <- vector(length = noPeaks)
            Afrags <- vector(length = noPeaks)
            for (j in 1:noPeaks) {
                Cfrags[j] = length(findFrag(UniqueUnitsOL$Frags[[i]]$Mass[j],
                  UniqueUnitsOL$Peaks)$Frag)
                Afrags[j] = length(findFrag(UniqueUnitsOL$Frags[[i]]$Mass[j],
                  allUnits$Peaks)$Frag)
                temp <- findFragMultiple((UniqueUnitsOL$Frags[[i]]$Mass[j] -
                  8), allUnits$Peaks, tol = 7.999)
                temp <- c(temp, findFragMultiple((UniqueUnitsOL$Frags[[i]]$Mass[j] +
                  8), allUnits$Peaks, tol = 7.999))
                if (length(temp) > 0)
                  for (k in 1:length(temp)) {
                    dist <- abs(UniqueUnitsOL$Frags[[i]]$Mass[j] -
                      temp[[k]]$Mass)
                    if (is.na(tClose[i]) || (tClose[i] > dist))
                      tClose[i] = dist
                  }
            }
        }
        fTable2 <- data.frame(tname, tBase, tDesc, tClose)
        names(fTable2) <- cnames
        write("\n\nInformative CpG_Units (Including Overlaps)",
            fname2, append = TRUE)
        write.csv(fTable2, fname2, append = TRUE)
        cat("Report written to:", fname2, "\n")
        options(warn = 0)
        tab <- peakTables(allUnits)
        for (i in nchar(fname):1) if (substr(fname, i, i) ==
            ".")
            break
        fname2 <- paste(substr(fname, 1, (i - 1)), " T Spectra.pdf",
            sep = "")
        pdf(fname2, width = 100, height = 20)
        op <- par("mai") * 2
        par(mai = op)
        plot(0, type = "n", xlim = c(minMass-100, maxMass+100), ylim = c(0,
            (max(tab[, 2] + tab[, 3]) + 5)), xlab = "Mass (Da)",
            ylab = "No Of Fragments", xaxp = c(minMass, maxMass, (maxMass-minMass)/100),
            yaxp = c(0, max(tab[, 2] + tab[, 3]), max(tab[, 2] +
                tab[, 3])), xaxs = "i", yaxs = "i", cex.axis = 2,
            cex.lab = 2, main = paste(fname, "T Cleavage Prediction"))
        for (i in 1:(nrow(tab) - 1)) {
            if ((tab[(i + 1), 1] - tab[i, 1]) < 15.9)
                rect(xleft = tab[i, 1], xright = tab[(i + 1),
                  1], ybottom = 0, ytop = (max(tab[, 2] + tab[,
                  3]) + 5), col = "yellow", border = NA)
        }
        rect(xleft = (tab[, 1] - 3), xright = (tab[, 1] + 3),
            ybottom = 0, ytop = tab[, 2], col = "green", border = NA)
        rect(xleft = (tab[, 1] - 3), xright = (tab[, 1] + 3),
            ybottom = tab[, 2], ytop = (tab[, 2] + tab[, 3]),
            col = "red", border = NA)
        for (i in 1:nrow(tab)) {
            lines(c(tab[i, 1], tab[i, 1]), c(sum(tab[i, 2:3]),
                sum(tab[i, 2:3], 0.5)), col = "purple")
            text((tab[i, 1] - 3), sum(tab[i, 2:3], 1), labels = tab[i,
                4], srt = 90, adj = 0, col = "green", cex = 0.8)
            text((tab[i, 1] + 3), sum(tab[i, 2:3], 1), labels = tab[i,
                5], srt = 90, adj = 0, col = "red", cex = 0.8)
        }
        legend("topright", c("Silent Fragments", "CpG Containing Fragments"),
            fill = c("green", "red"), )
        dev.off()
        cat("Predicted spectra written to:", fname2, "\n")
        ampSize <- AMPa[[length(AMPa)]]$Pos+AMPa[[length(AMPa)]]$Len
        CpGPositions <- vector()
        CpGCol <- vector()
        for (i in 1:length(CpGUnits$Frags)) {
          for (j in nchar(CpGUnits$Frags[[i]]$Seq):1) {
            if (substr(CpGUnits$Frags[[i]]$Seq,j,j)=="*") {
              CpGPositions = c(CpGPositions, ampSize-(CpGUnits$Frags[[i]]$Pos+j))
              if (CpGUnits$Frags[[i]]$Mass[1]>minMass&&CpGUnits$Frags[[i]]$Mass[length(CpGUnits$Frags[[i]]$Mass)]<maxMass) CpGCol = c(CpGCol,"red") else CpGCol = c(CpGCol, "gray")
              }
          }
        }
        cleaves <- vector()
        cleavesCol <- vector()
        for (i in 1:length(AMPa)) {
          cleaves = c(cleaves,(ampSize-(AMPa[[i]]$Pos+AMPa[[i]]$Len)))
          if (AMPa[[i]]$Mass>minMass&&AMPa[[i]]$Mass<maxMass) cleavesCol = c(cleavesCol, "black") else cleavesCol = c(cleavesCol, "gray")
        }
        for (i in nchar(fname):1) if (substr(fname, i, i) ==
            ".")
            break
        fname3 <- paste(substr(fname, 1, (i - 1)), " T Fragmentation.pdf",
            sep = "")        
        pdf(fname3, width=20, height=10, version="1.4")
        plot(0,xlim=c(0,ampSize),ylim=c(-1.1,1.1),ann=FALSE, xaxt="n", yaxt="n", type="n")
        segments(c(ampSize,cleaves), 0.8, c(cleaves,0), 0.8, c(cleavesCol, "black"))
        segments(c(10,ampSize-11), c(0.75,0.75),c(10,ampSize-11), c(0.9,0.9))
        text(x=c(10,ampSize-11), y=(0.95), labels=paste(c(1,(ampSize-22)),"bp"), cex=1.5)
        points(CpGPositions, rep(0.8,length(CpGPositions)), pch=19, col=CpGCol, cex=0.9)
        text(x=CpGPositions,y=0.9,labels=c(paste("CpG",1:length(CpGPositions))), cex=0.9, srt=90, adj=0)
        segments((cleaves+0.01),0.78,(cleaves-0.01),0.82)
        rect(xleft=0, ybottom=0.78, xright=10, ytop=0.82, col="#FFFF0099", border=NA)
        rect(xleft=ampSize-11, ybottom=0.78, xright=ampSize, ytop=0.82, col="#FFFF0099", border=NA)
        unitsPos <- vector()
        for (i in 1:length(CpGUnits$Frags)) {
          CpGPos <- vector()
          for (j in 1:nchar(CpGUnits$Frags[[i]]$Seq)) {
            if (substr(CpGUnits$Frags[[i]]$Seq,j,j)=="*")
              CpGPos = c(CpGPos, ampSize-(CpGUnits$Frags[[i]]$Pos+j))
          }
          segments(CpGPos, 0.77, mean(CpGPos), 0.7, col="violet")
          segments(mean(CpGPos), 0.7, mean(CpGPos), 0.6, col="violet")
          text(x=mean(CpGPos), y=0.58, label=CpGUnits$Frags[[i]]$Name, srt=90, cex=1.0, adj=1)
          unitsPos = c(unitsPos, mean(CpGPos))
        }
        OLs <- vector()
        OLfrom <- vector()
        for (j in 1:length(UniqueUnitsOL$Frags))  {
          if (length(UniqueUnitsOL$Frags[[j]]$FromFrag)>1) {
            OLs <- c(OLs,UniqueUnitsOL$Frags[[j]]$FromFrag[2:length(UniqueUnitsOL$Frags[[j]]$FromFrag)])
            OLfrom <- c(OLfrom, rep(UniqueUnitsOL$Frags[[j]]$FromFrag[1],(length(UniqueUnitsOL$Frags[[j]]$FromFrag)-1))) 
          }
        }
        OLcol <- rainbow(n=length(OLs))
        for (i in 1:length(UniqueUnitsOL$Frags)) {
          segments(unitsPos[UniqueUnitsOL$Frags[[i]]$FromFrag[1]],0.25,unitsPos[UniqueUnitsOL$Frags[[i]]$FromFrag[1]],0.2-(0.05*length(OLs)),col="violet")
          text(x=unitsPos[UniqueUnitsOL$Frags[[i]]$FromFrag[1]], y=0.18-(0.05*length(OLs)), label=paste(UniqueUnitsOL$Frags[[i]]$Name," (",UniqueUnitsOL$Frags[[i]]$CpGSites,")",sep=""), srt=90, cex=1, adj=1)
        }
        if (length(OLs)>0) for (i in 1:length(OLs)) {
          segments(unitsPos[OLs[i]],0.25,unitsPos[OLs[i]],0.2-(0.05*(i-1)), col=OLcol[i])
          segments(unitsPos[OLs[i]],0.2-(0.05*(i-1)),unitsPos[OLfrom[i]],0.2-(0.05*(i-1)), col=OLcol[i])
        }
        text(x=ampSize/2, y=-1.1, label=paste(fname, "T Cleavage reaction predicted fragmentation"), adj=0.5, cex=1.5)
        dev.off()
        cat("Predicted fragmentation written to:", fname3, "\n")
    }
    showCleavage <- function(amp) {
        temp <- vector()
        for (i in 1:length(amp)) temp = paste(temp, amp[[i]]$Seq,
            sep = ",")
        for (i in 0:(nchar(temp) - 1)) {
            t2 <- substr(temp, (nchar(temp) - i), (nchar(temp) -
                i))
            if (i == 0) {
                bottom <- t2
            }
            else bottom <- paste(bottom, t2, sep = "")
        }
        top <- gsub("g", "C", gsub("t", "A", gsub("a", "T", gsub("c",
            "G", bottom))))
        top <- gsub("C", "c", gsub("T", "t", gsub("A", "a", gsub("G",
            "g", top))))
        top = paste("5'-", top, "-3'", sep = ",")
        bottom = paste("3'-", bottom, "-5'", sep = ",")
        return(c(top, bottom))
    }
    makeCleave <- function(seq, pos) {
        CpGs <- 0
        Amass <- 329.2098
        Gmass <- 345.2091
        Cmass <- 289.1851
        Tmass <- 306.169
        mass <- 19.02327
        if (pos == 1)
            mass <- mass + 239.939
        for (j in 1:nchar(seq)) {
            if (substr(seq, j, j) == "a") {
                mass[[1]] = mass[[1]] + Amass
            }
            if (substr(seq, j, j) == "c") {
                mass[[1]] = mass[[1]] + Cmass
            }
            if (substr(seq, j, j) == "t") {
                mass[[1]] = mass[[1]] + Tmass
            }
            if (substr(seq, j, j) == "*") {
                mass[[1]] = mass[[1]] + Amass
                CpGs = CpGs + 1
            }
            if (substr(seq, j, j) == "g") {
                mass[[1]] = mass[[1]] + Gmass
            }
        }
        if (CpGs != 0)
            for (j in 1:CpGs) {
                mass = c(mass, (mass[[length(mass)]] + (Gmass -
                  Amass)))
            }
        return(list(Seq = seq, Pos = pos, Len = nchar(seq), CpGs = CpGs,
            Mass = round(mass, 3)))
    }
    parseSequence <- function(filename) {
        seq <- paste(toupper(scan(filename, what = "character",
            sep = "\n", quiet = TRUE)))
        if (length(seq) > 1) {
            seq2 <- seq
            seq <- seq[1]
            for (i in 2:length(seq2)) seq = paste(seq[1], seq2[i])
        }
        seq <- gsub(" ", "", seq)
        cpg <- gsub("CG", "*G", seq)
        bis <- paste("AGGAAGAGAG", gsub("C", "T", cpg), "AGCCTTCTCCC",
            sep = "")
        temp <- gsub("G", "c", gsub("T", "a", gsub("A", "t",
            gsub("C", "g", bis))))
        for (i in 0:(nchar(temp) - 1)) {
            t2 <- substr(temp, (nchar(temp) - i), (nchar(temp) -
                i))
            if (i == 0) {
                revcom <- t2
            }
            else revcom <- paste(revcom, t2, sep = "")
        }
        upto <- 0
        for (i in 1:nchar(revcom)) {
            if (substr(revcom, i, i) == "t") {
                if (upto == 0) {
                  cleaves <- list(makeCleave(substr(revcom, 1,
                    i), 1))
                }
                else cleaves <- c(cleaves, list(makeCleave(substr(revcom,
                  (upto + 1), i), (upto + 1))))
                upto = i
            }
        }
        if (upto != (nchar(revcom)))
            cleaves <- c(cleaves, list(makeCleave(substr(revcom,
                (upto + 1), nchar(revcom)), (upto + 1))))
        return(cleaves)
    }
    makePeaks <- function(frags) {
        first <- TRUE
        for (i in 1:length(frags)) for (j in 1:length(frags[[i]]$Mass)) {
            if (first == TRUE) {
                first = FALSE
                peaks = list(list(Mass = frags[[i]]$Mass[j],
                  Frag = i, Peak = j))
            }
            else {
                dup <- 0
                for (k in 1:length(peaks)) if (identical(peaks[[k]]$Mass,
                  frags[[i]]$Mass[j]))
                  dup = k
                if (dup == 0) {
                  peaks <- c(peaks, list(list(Mass = frags[[i]]$Mass[j],
                    Frag = i, Peak = j)))
                }
                else {
                  peaks[[dup]]$Frag = c(peaks[[dup]]$Frag, i)
                  peaks[[dup]]$Peak = c(peaks[[dup]]$Peak, j)
                }
            }
        }
        return(peaks)
    }
    findPeak <- function(pMass, pTable, tol = 1) {
        for (i in 1:length(pTable$Reference.mass)) if ((pTable$Reference.mass[i] +
            tol) > pMass && (pTable$Reference.mass[i] - tol) <
            pMass)
            return(pTable[i, ])
        return(list(SNR = NA))
    }
    findUniqueUnits <- function(cleaves, minMass, maxMass,
        skipOL = FALSE) {
        units <- list()
        for (i in 1:length(cleaves)) if (cleaves[[i]]$Mass[1] >
            minMass && cleaves[[i]]$Mass[length(cleaves[[i]]$Mass)] <
            maxMass)
            if (length(units) == 0) {
                units <- cleaves[i]
                units[[1]] = c(units[[1]], list(FromFrag = i))
            }
            else {
                dup <- 0
                for (j in 1:length(units)) if (cleaves[[i]]$Mass[1] ==
                  units[[j]]$Mass[1])
                  dup <- j
                if (dup == 0) {
                  units <- c(units, cleaves[i])
                  units[[length(units)]] = c(units[[length(units)]],
                    list(FromFrag = i))
                }
                else if (length(cleaves[[i]]$Mass) > length(units[[dup]]$Mass)) {
                  temp <- units[[dup]]
                  units[[dup]] = cleaves[[i]]
                  units[[dup]]$FromFrag = c(i,temp$FromFrag)
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites,
                    "_OL_", temp$CpGSites, sep = "")
                }
                else if (length(cleaves[[i]]$Mass) < length(units[[dup]]$Mass)) {
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites,
                    "_OL_", cleaves[[i]]$CpGSites, sep = "")
                  units[[dup]]$FromFrag = c(units[[dup]]$FromFrag,i)                    
                }
                else {
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites,
                    "_D_", cleaves[[i]]$CpGSites, sep = "")
                  units[[dup]]$FromFrag = c(units[[dup]]$FromFrag,i)                    
                }
            }
        if (skipOL) {
            unitsNew <- list()
            for (i in 1:length(units)) if (length(grep("_OL_",
                units[[i]]$CpGSites)) == 0)
                unitsNew = c(unitsNew, units[i])
            units = unitsNew
        }
        return(units)
    }
    findFrag <- function(pMass, peaks, tol = 1) {
        for (i in 1:length(peaks)) if ((peaks[[i]]$Mass + tol) >
            pMass && (peaks[[i]]$Mass - tol) < pMass)
            return(peaks[[i]])
        return(list(Frag = NA))
    }
    findFragMultiple <- function(pMass, peaks, tol = 1) {
        frags <- list()
        for (i in 1:length(peaks)) if ((peaks[[i]]$Mass + tol) >
            pMass && (peaks[[i]]$Mass - tol) < pMass)
            frags = c(frags, peaks[i])
        return(frags)
    }
    findCpGs <- function(cleaves) {
        numUnits <- 0
        noCpGs <- 0
        for (i in length(cleaves):1) {
            if (cleaves[[i]]$CpGs > 0) {
                if (numUnits == 0)
                  units <- cleaves[i]
                else units <- c(units, cleaves[i])
                numUnits = numUnits + 1
                CpGName <- "CpG"
                for (j in (noCpGs + 1):(noCpGs + cleaves[[i]]$CpGs)) CpGName = paste(CpGName,
                  j, sep = "_")
                units[[numUnits]] = c(units[[numUnits]], list(Name = paste("CpG_Unit_",
                  numUnits, sep = ""), CpGSites = CpGName))
                noCpGs = noCpGs + cleaves[[i]]$CpGs
            }
        }
        return(units)
    }
    "peakTables" <- function(units) {
        noPeaks <- length(units$Peaks)
        mass <- vector(mode = "numeric")
        silent <- vector(mode = "numeric")
        meth <- vector(mode = "numeric")
        silentD <- vector(mode = "character")
        methD <- vector(mode = "character")
        for (i in 1:noPeaks) if (units$Peaks[[i]]$Mass > minMass &&
            units$Peaks[[i]]$Mass < maxMass) {
            mass = c(mass, units$Peaks[[i]]$Mass)
            meth = c(meth, 0)
            silent = c(silent, 0)
            methD = c(methD, "")
            silentD = c(silentD, "")
            ii <- length(mass)
            for (j in 1:length(units$Peaks[[i]]$Frag)) if (units$Frags[[units$Peaks[[i]]$Frag[j]]]$CpG >
                0) {
                meth[ii] = meth[ii] + 1
                tDesc <- paste(units$Frags[[units$Peaks[[i]]$Frag[j]]]$Seq,
                  "@", units$Frags[[units$Peaks[[i]]$Frag[j]]]$Pos,
                  sep = "")
                methD[ii] = paste(methD[ii], tDesc)
            }
            else {
                silent[ii] = silent[ii] + 1
                tDesc <- paste(units$Frags[[units$Peaks[[i]]$Frag[j]]]$Seq,
                  "@", units$Frags[[units$Peaks[[i]]$Frag[j]]]$Pos,
                  sep = "")
                silentD[ii] = paste(silentD[ii], tDesc)
            }
        }
        o <- order(mass)
        table <- data.frame(mass[o], silent[o], meth[o], silentD[o],
            methD[o])
        return(table)
    }
    if (length(fnames) > 0)
        for (i in 1:length(fnames)) generateReport(fnames[i])
}

ampliconReportC <- function (fnames = NA, minMass, maxMass)
{
    generateReport <- function(fname) {
        AMPa <- parseSequence(fname)
        AMPc <- findCpGs(AMPa)
        AMPu <- findUniqueUnits(AMPc, minMass, maxMass, skipOL = TRUE)
        AMPuOL <- findUniqueUnits(AMPc, minMass, maxMass, skipOL = FALSE)
        allUnits <- list(Frags = AMPa, Peaks = makePeaks(AMPa))
        CpGUnits <- list(Frags = AMPc, Peaks = makePeaks(AMPc))
        UniqueUnits <- list(Frags = AMPu, Peaks = makePeaks(AMPc))
        UniqueUnitsOL <- list(Frags = AMPuOL, Peaks = makePeaks(AMPc))
        header <- rep("", 11)
        header[1] = paste("C Cleavage report of amplicon", fname)
        header[2] = paste("No of Fragments:", length(allUnits$Frags))
        header[3] = paste("No of Peaks:", length(allUnits$Peaks))
        CpGCount <- 0
        nCpGUnits <- length(CpGUnits$Frags)
        for (i in 1:length(CpGUnits$Frags)) CpGCount = CpGCount +
            CpGUnits$Frags[[i]]$CpGs
        header[4] = paste(CpGCount, "CpGs", "in", nCpGUnits,
            "CpG Units")
        CpGCount <- 0
        nUniqueUnits <- length(UniqueUnitsOL$Frags)
        for (i in 1:length(UniqueUnitsOL$Frags)) CpGCount = CpGCount +
            UniqueUnitsOL$Frags[[i]]$CpGs
        header[5] = paste(CpGCount, "CpGs", "in", nUniqueUnits,
            "informative CpG Units (including overlaps)")
        header[7] = "Predicted C Cleavage fragments - top (template) and bottom (transcribed) strand"
        header[8:9] = showCleavage(AMPa)
        for (i in nchar(fname):1) if (substr(fname, i, i) ==
            ".")
            break
        fname2 <- paste(substr(fname, 1, (i - 1)), " C Report.csv",
            sep = "")
        options(warn = -1)
        write(header, fname2, append = FALSE)
        nfrags <- length(CpGUnits$Frags)
        cnames <- c("Name", "Position", "Sequence", "CpGs", "Base mass",
            "Desc", "<15.9Da")
        tname <- vector(mode = "character", length = nfrags)
        tPos <- vector(mode = "numeric", length = nfrags)
        tSeq <- vector(mode = "character", length = nfrags)
        tNo <- vector(mode = "numeric", length = nfrags)
        tBase <- vector(mode = "numeric", length = nfrags)
        tDesc <- vector(mode = "character", length = nfrags)
        tClose <- rep(NA, nfrags)
        for (i in 1:nfrags) {
            tname[i] = CpGUnits$Frags[[i]]$Name
            tNo[i] = CpGUnits$Frags[[i]]$CpGs
            tPos[i] = CpGUnits$Frags[[i]]$Pos
            tSeq[i] = CpGUnits$Frags[[i]]$Seq
            tBase[i] = CpGUnits$Frags[[i]]$Mass[1]
            tDesc[i] = CpGUnits$Frags[[i]]$CpGSites
            noPeaks <- length(CpGUnits$Frags[[i]]$Mass)
            Cfrags <- vector(length = noPeaks)
            Afrags <- vector(length = noPeaks)
            for (j in 1:noPeaks) {
                Cfrags[j] = length(findFrag(CpGUnits$Frags[[i]]$Mass[j],
                  CpGUnits$Peaks)$Frag)
                Afrags[j] = length(findFrag(CpGUnits$Frags[[i]]$Mass[j],
                  allUnits$Peaks)$Frag)
                temp <- findFragMultiple((CpGUnits$Frags[[i]]$Mass[j] -
                  8), allUnits$Peaks, tol = 7.999)
                temp <- c(temp, findFragMultiple((CpGUnits$Frags[[i]]$Mass[j] +
                  8), allUnits$Peaks, tol = 7.999))
                if (length(temp) > 0)
                  for (k in 1:length(temp)) {
                    dist <- abs(CpGUnits$Frags[[i]]$Mass[j] -
                      temp[[k]]$Mass)
                    if (is.na(tClose[i]) || (tClose[i] > dist))
                      tClose[i] = dist
                  }
            }
        }
        fTable <- data.frame(tname, tPos, tSeq, tNo, tBase, tDesc,
            tClose)
        names(fTable) <- cnames
        write("All CpG_Units", fname2, append = TRUE)
        write.csv(fTable, fname2, append = TRUE)
        nfrags <- length(UniqueUnitsOL$Frags)
        cnames <- c("Name", "Base mass", "Desc", "<15.9Da")
        tname <- vector(mode = "character", length = nfrags)
        tBase <- vector(mode = "numeric", length = nfrags)
        tDesc <- vector(mode = "character", length = nfrags)
        tClose <- rep(NA, nfrags)
        for (i in 1:nfrags) {
            tname[i] = UniqueUnitsOL$Frags[[i]]$Name
            tBase[i] = UniqueUnitsOL$Frags[[i]]$Mass[1]
            tDesc[i] = UniqueUnitsOL$Frags[[i]]$CpGSites
            noPeaks <- length(UniqueUnitsOL$Frags[[i]]$Mass)
            Cfrags <- vector(length = noPeaks)
            Afrags <- vector(length = noPeaks)
            for (j in 1:noPeaks) {
                Cfrags[j] = length(findFrag(UniqueUnitsOL$Frags[[i]]$Mass[j],
                  UniqueUnitsOL$Peaks)$Frag)
                Afrags[j] = length(findFrag(UniqueUnitsOL$Frags[[i]]$Mass[j],
                  allUnits$Peaks)$Frag)
                temp <- findFragMultiple((UniqueUnitsOL$Frags[[i]]$Mass[j] -
                  8), allUnits$Peaks, tol = 7.999)
                temp <- c(temp, findFragMultiple((UniqueUnitsOL$Frags[[i]]$Mass[j] +
                  8), allUnits$Peaks, tol = 7.999))
                if (length(temp) > 0)
                  for (k in 1:length(temp)) {
                    dist <- abs(UniqueUnitsOL$Frags[[i]]$Mass[j] -
                      temp[[k]]$Mass)
                    if (is.na(tClose[i]) || (tClose[i] > dist))
                      tClose[i] = dist
                  }
            }
        }
        fTable2 <- data.frame(tname, tBase, tDesc, tClose)
        names(fTable2) <- cnames
        write("\n\nInformative CpG_Units (Including Overlaps)",
            fname2, append = TRUE)
        write.csv(fTable2, fname2, append = TRUE)
        cat("Report written to:", fname2, "\n")
        options(warn = 0)
        tab <- peakTables(allUnits)
        for (i in nchar(fname):1) if (substr(fname, i, i) ==
            ".")
            break
        fname2 <- paste(substr(fname, 1, (i - 1)), " C Spectra.pdf",
            sep = "")
        pdf(fname2, width = 100, height = 20)
        op <- par("mai") * 2
        par(mai = op)
        plot(0, type = "n", xlim = c(minMass-100, maxMass+100), ylim = c(0,
            (max(tab[, 2] + tab[, 3]) + 5)), xlab = "Mass (Da)",
            ylab = "No Of Fragments", xaxp = c(minMass, maxMass, (maxMass-minMass)/100),
            yaxp = c(0, max(tab[, 2] + tab[, 3]), max(tab[, 2] +
                tab[, 3])), xaxs = "i", yaxs = "i", cex.axis = 2,
            cex.lab = 2, main = paste(fname, "T Cleavage Prediction"))
        for (i in 1:(nrow(tab) - 1)) {
            if ((tab[(i + 1), 1] - tab[i, 1]) < 15.9)
                rect(xleft = tab[i, 1], xright = tab[(i + 1),
                  1], ybottom = 0, ytop = (max(tab[, 2] + tab[,
                  3]) + 5), col = "yellow", border = NA)
        }
        rect(xleft = (tab[, 1] - 3), xright = (tab[, 1] + 3),
            ybottom = 0, ytop = tab[, 2], col = "green", border = NA)
        rect(xleft = (tab[, 1] - 3), xright = (tab[, 1] + 3),
            ybottom = tab[, 2], ytop = (tab[, 2] + tab[, 3]),
            col = "red", border = NA)
        for (i in 1:nrow(tab)) {
            lines(c(tab[i, 1], tab[i, 1]), c(sum(tab[i, 2:3]),
                sum(tab[i, 2:3], 0.5)), col = "purple")
            text((tab[i, 1] - 3), sum(tab[i, 2:3], 1), labels = tab[i,
                4], srt = 90, adj = 0, col = "green", cex = 0.8)
            text((tab[i, 1] + 3), sum(tab[i, 2:3], 1), labels = tab[i,
                5], srt = 90, adj = 0, col = "red", cex = 0.8)
        }
        legend("topright", c("Silent Fragments", "CpG Containing Fragments"),
            fill = c("green", "red"), )
        dev.off()
        cat("Predicted spectra written to:", fname2, "\n")
        ampSize <- AMPa[[length(AMPa)]]$Pos+AMPa[[length(AMPa)]]$Len
        CpGPositions <- vector()
        CpGCol <- vector()
        for (i in 1:length(CpGUnits$Frags)) {
          for (j in nchar(CpGUnits$Frags[[i]]$Seq):1) {
            if (substr(CpGUnits$Frags[[i]]$Seq,j,j)=="*") {
              CpGPositions = c(CpGPositions, ampSize-(CpGUnits$Frags[[i]]$Pos+j))
              if (CpGUnits$Frags[[i]]$Mass[1]>minMass&&CpGUnits$Frags[[i]]$Mass[length(CpGUnits$Frags[[i]]$Mass)]<maxMass) CpGCol = c(CpGCol,"red") else CpGCol = c(CpGCol, "gray")
              }
          }
        }
        cleaves <- vector()
        cleavesCol <- vector()
        for (i in 1:length(AMPa)) {
          cleaves = c(cleaves,(ampSize-(AMPa[[i]]$Pos+AMPa[[i]]$Len)))
          if (AMPa[[i]]$Mass>minMass&&AMPa[[i]]$Mass<maxMass) cleavesCol = c(cleavesCol, "black") else cleavesCol = c(cleavesCol, "gray")
        }
        for (i in nchar(fname):1) if (substr(fname, i, i) ==
            ".")
            break
        fname3 <- paste(substr(fname, 1, (i - 1)), " C Fragmentation.pdf",
            sep = "")        
        pdf(fname3, width=20, height=10, version="1.4")
        plot(0,xlim=c(0,ampSize),ylim=c(-1.1,1.1),ann=FALSE, xaxt="n", yaxt="n", type="n")
        segments(c(ampSize,cleaves), 0.8, c(cleaves,0), 0.8, c(cleavesCol, "black"))
        segments(c(10,ampSize-11), c(0.75,0.75),c(10,ampSize-11), c(0.9,0.9))
        text(x=c(10,ampSize-11), y=(0.95), labels=paste(c(1,(ampSize-22)),"bp"), cex=1.5)
        points(CpGPositions, rep(0.8,length(CpGPositions)), pch=19, col=CpGCol, cex=0.9)
        text(x=CpGPositions,y=0.9,labels=c(paste("CpG",1:length(CpGPositions))), cex=0.9, srt=90, adj=0)
        segments((cleaves+0.01),0.78,(cleaves-0.01),0.82)
        rect(xleft=0, ybottom=0.78, xright=10, ytop=0.82, col="#FFFF0099", border=NA)
        rect(xleft=ampSize-11, ybottom=0.78, xright=ampSize, ytop=0.82, col="#FFFF0099", border=NA)
        unitsPos <- vector()
        for (i in 1:length(CpGUnits$Frags)) {
          CpGPos <- vector()
          for (j in 1:nchar(CpGUnits$Frags[[i]]$Seq)) {
            if (substr(CpGUnits$Frags[[i]]$Seq,j,j)=="*")
              CpGPos = c(CpGPos, ampSize-(CpGUnits$Frags[[i]]$Pos+j))
          }
          segments(CpGPos, 0.77, mean(CpGPos), 0.7, col="violet")
          segments(mean(CpGPos), 0.7, mean(CpGPos), 0.6, col="violet")
          text(x=mean(CpGPos), y=0.58, label=CpGUnits$Frags[[i]]$Name, srt=90, cex=1.0, adj=1)
          unitsPos = c(unitsPos, mean(CpGPos))
        }
        OLs <- vector()
        OLfrom <- vector()
        for (j in 1:length(UniqueUnitsOL$Frags))  {
          if (length(UniqueUnitsOL$Frags[[j]]$FromFrag)>1) {
            OLs <- c(OLs,UniqueUnitsOL$Frags[[j]]$FromFrag[2:length(UniqueUnitsOL$Frags[[j]]$FromFrag)])
            OLfrom <- c(OLfrom, rep(UniqueUnitsOL$Frags[[j]]$FromFrag[1],(length(UniqueUnitsOL$Frags[[j]]$FromFrag)-1))) 
          }
        }
        OLcol <- rainbow(n=length(OLs))
        for (i in 1:length(UniqueUnitsOL$Frags)) {
          segments(unitsPos[UniqueUnitsOL$Frags[[i]]$FromFrag[1]],0.25,unitsPos[UniqueUnitsOL$Frags[[i]]$FromFrag[1]],0.2-(0.05*length(OLs)),col="violet")
          text(x=unitsPos[UniqueUnitsOL$Frags[[i]]$FromFrag[1]], y=0.18-(0.05*length(OLs)), label=paste(UniqueUnitsOL$Frags[[i]]$Name," (",UniqueUnitsOL$Frags[[i]]$CpGSites,")",sep=""), srt=90, cex=1, adj=1)
        }
        if (length(OLs)>0) for (i in 1:length(OLs)) {
          segments(unitsPos[OLs[i]],0.25,unitsPos[OLs[i]],0.2-(0.05*(i-1)), col=OLcol[i])
          segments(unitsPos[OLs[i]],0.2-(0.05*(i-1)),unitsPos[OLfrom[i]],0.2-(0.05*(i-1)), col=OLcol[i])
        }
        text(x=ampSize/2, y=-1.1, label=paste(fname, "C Cleavage reaction predicted fragmentation"), adj=0.5, cex=1.5)
        dev.off()
        cat("Predicted fragmentation written to:", fname3, "\n")
    }
    showCleavage <- function(amp) {
        temp <- vector()
        for (i in 1:length(amp)) temp = paste(temp, amp[[i]]$Seq,
            sep = ",")
        for (i in 0:(nchar(temp) - 1)) {
            t2 <- substr(temp, (nchar(temp) - i), (nchar(temp) -
                i))
            if (i == 0) {
                bottom <- t2
            }
            else bottom <- paste(bottom, t2, sep = "")
        }
        top <- gsub("g", "C", gsub("t", "A", gsub("a", "T", gsub("c",
            "G", bottom))))
        top <- gsub("C", "c", gsub("T", "t", gsub("A", "a", gsub("G",
            "g", top))))
        top = paste("5'-", top, "-3'", sep = ",")
        bottom = paste("3'-", bottom, "-5'", sep = ",")
        return(c(top, bottom))
    }
    makeCleave <- function(seq, pos) {
        CpGs <- 0
        Amass <- 329.2098
        Gmass <- 345.2091
        Cmass <- 305.1844
        Tmass <- 304.1967
        mass <- 19.02327
        if (pos == 1)
            mass <- mass + 239.939
        for (j in 1:nchar(seq)) {
            if (substr(seq, j, j) == "a") {
                mass[[1]] = mass[[1]] + Amass
            }
            if (substr(seq, j, j) == "c") {
                mass[[1]] = mass[[1]] + Cmass
            }
            if (substr(seq, j, j) == "t") {
                mass[[1]] = mass[[1]] + Tmass
            }
            if (substr(seq, j, j) == "*") {
                mass[[1]] = mass[[1]] + Amass
                CpGs = CpGs + 1
            }
            if (substr(seq, j, j) == "g") {
                mass[[1]] = mass[[1]] + Gmass
            }
        }
        if (CpGs != 0)
            for (j in 1:CpGs) {
                mass = c(mass, (mass[[length(mass)]] + (Gmass -
                  Amass)))
            }
        return(list(Seq = seq, Pos = pos, Len = nchar(seq), CpGs = CpGs,
            Mass = round(mass, 3)))
    }
    parseSequence <- function(filename) {
        seq <- paste(toupper(scan(filename, what = "character",
            sep = "\n", quiet = TRUE)))
        if (length(seq) > 1) {
            seq2 <- seq
            seq <- seq[1]
            for (i in 2:length(seq2)) seq = paste(seq[1], seq2[i])
        }
        seq <- gsub(" ", "", seq)
        cpg <- gsub("CG", "*G", seq)
        bis <- paste("AGGAAGAGAG", gsub("C", "T", cpg), "AGCCTTCTCCC",
            sep = "")
        temp <- gsub("G", "c", gsub("T", "a", gsub("A", "t",
            gsub("C", "g", bis))))
        for (i in 0:(nchar(temp) - 1)) {
            t2 <- substr(temp, (nchar(temp) - i), (nchar(temp) -
                i))
            if (i == 0) {
                revcom <- t2
            }
            else revcom <- paste(revcom, t2, sep = "")
        }
        upto <- 0
        for (i in 1:nchar(revcom)) {
            if (substr(revcom, i, i) == "c") {
                if (upto == 0) {
                  cleaves <- list(makeCleave(substr(revcom, 1,
                    i), 1))
                }
                else cleaves <- c(cleaves, list(makeCleave(substr(revcom,
                  (upto + 1), i), (upto + 1))))
                upto = i
            }
        }
        if (upto != (nchar(revcom)))
            cleaves <- c(cleaves, list(makeCleave(substr(revcom,
                (upto + 1), nchar(revcom)), (upto + 1))))
        return(cleaves)
    }
    makePeaks <- function(frags) {
        first <- TRUE
        for (i in 1:length(frags)) for (j in 1:length(frags[[i]]$Mass)) {
            if (first == TRUE) {
                first = FALSE
                peaks = list(list(Mass = frags[[i]]$Mass[j],
                  Frag = i, Peak = j))
            }
            else {
                dup <- 0
                for (k in 1:length(peaks)) if (identical(peaks[[k]]$Mass,
                  frags[[i]]$Mass[j]))
                  dup = k
                if (dup == 0) {
                  peaks <- c(peaks, list(list(Mass = frags[[i]]$Mass[j],
                    Frag = i, Peak = j)))
                }
                else {
                  peaks[[dup]]$Frag = c(peaks[[dup]]$Frag, i)
                  peaks[[dup]]$Peak = c(peaks[[dup]]$Peak, j)
                }
            }
        }
        return(peaks)
    }
    findPeak <- function(pMass, pTable, tol = 1) {
        for (i in 1:length(pTable$Reference.mass)) if ((pTable$Reference.mass[i] +
            tol) > pMass && (pTable$Reference.mass[i] - tol) <
            pMass)
            return(pTable[i, ])
        return(list(SNR = NA))
    }
    findUniqueUnits <- function(cleaves, minMass, maxMass,
        skipOL = FALSE) {
        units <- list()
        for (i in 1:length(cleaves)) if (cleaves[[i]]$Mass[1] >
            minMass && cleaves[[i]]$Mass[length(cleaves[[i]]$Mass)] <
            maxMass)
            if (length(units) == 0) {
                units <- cleaves[i]
                units[[1]] = c(units[[1]], list(FromFrag = i))
            }
            else {
                dup <- 0
                for (j in 1:length(units)) if (cleaves[[i]]$Mass[1] ==
                  units[[j]]$Mass[1])
                  dup <- j
                if (dup == 0) {
                  units <- c(units, cleaves[i])
                  units[[length(units)]] = c(units[[length(units)]],
                    list(FromFrag = i))
                }
                else if (length(cleaves[[i]]$Mass) > length(units[[dup]]$Mass)) {
                  temp <- units[[dup]]
                  units[[dup]] = cleaves[[i]]
                  units[[dup]]$FromFrag = c(i,temp$FromFrag)
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites,
                    "_OL_", temp$CpGSites, sep = "")
                }
                else if (length(cleaves[[i]]$Mass) < length(units[[dup]]$Mass)) {
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites,
                    "_OL_", cleaves[[i]]$CpGSites, sep = "")
                  units[[dup]]$FromFrag = c(units[[dup]]$FromFrag,i)                    
                }
                else {
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites,
                    "_D_", cleaves[[i]]$CpGSites, sep = "")
                  units[[dup]]$FromFrag = c(units[[dup]]$FromFrag,i)                    
                }
            }
        if (skipOL) {
            unitsNew <- list()
            for (i in 1:length(units)) if (length(grep("_OL_",
                units[[i]]$CpGSites)) == 0)
                unitsNew = c(unitsNew, units[i])
            units = unitsNew
        }
        return(units)
    }
    findFrag <- function(pMass, peaks, tol = 1) {
        for (i in 1:length(peaks)) if ((peaks[[i]]$Mass + tol) >
            pMass && (peaks[[i]]$Mass - tol) < pMass)
            return(peaks[[i]])
        return(list(Frag = NA))
    }
    findFragMultiple <- function(pMass, peaks, tol = 1) {
        frags <- list()
        for (i in 1:length(peaks)) if ((peaks[[i]]$Mass + tol) >
            pMass && (peaks[[i]]$Mass - tol) < pMass)
            frags = c(frags, peaks[i])
        return(frags)
    }
    findCpGs <- function(cleaves) {
        numUnits <- 0
        noCpGs <- 0
        for (i in length(cleaves):1) {
            if (cleaves[[i]]$CpGs > 0) {
                if (numUnits == 0)
                  units <- cleaves[i]
                else units <- c(units, cleaves[i])
                numUnits = numUnits + 1
                CpGName <- "CpG"
                for (j in (noCpGs + 1):(noCpGs + cleaves[[i]]$CpGs)) CpGName = paste(CpGName,
                  j, sep = "_")
                units[[numUnits]] = c(units[[numUnits]], list(Name = paste("CpG_Unit_",
                  numUnits, sep = ""), CpGSites = CpGName))
                noCpGs = noCpGs + cleaves[[i]]$CpGs
            }
        }
        return(units)
    }
    "peakTables" <- function(units) {
        noPeaks <- length(units$Peaks)
        mass <- vector(mode = "numeric")
        silent <- vector(mode = "numeric")
        meth <- vector(mode = "numeric")
        silentD <- vector(mode = "character")
        methD <- vector(mode = "character")
        for (i in 1:noPeaks) if (units$Peaks[[i]]$Mass > minMass &&
            units$Peaks[[i]]$Mass < maxMass) {
            mass = c(mass, units$Peaks[[i]]$Mass)
            meth = c(meth, 0)
            silent = c(silent, 0)
            methD = c(methD, "")
            silentD = c(silentD, "")
            ii <- length(mass)
            for (j in 1:length(units$Peaks[[i]]$Frag)) if (units$Frags[[units$Peaks[[i]]$Frag[j]]]$CpG >
                0) {
                meth[ii] = meth[ii] + 1
                tDesc <- paste(units$Frags[[units$Peaks[[i]]$Frag[j]]]$Seq,
                  "@", units$Frags[[units$Peaks[[i]]$Frag[j]]]$Pos,
                  sep = "")
                methD[ii] = paste(methD[ii], tDesc)
            }
            else {
                silent[ii] = silent[ii] + 1
                tDesc <- paste(units$Frags[[units$Peaks[[i]]$Frag[j]]]$Seq,
                  "@", units$Frags[[units$Peaks[[i]]$Frag[j]]]$Pos,
                  sep = "")
                silentD[ii] = paste(silentD[ii], tDesc)
            }
        }
        o <- order(mass)
        table <- data.frame(mass[o], silent[o], meth[o], silentD[o],
            methD[o])
        return(table)
    }
    if (length(fnames) > 0)
        for (i in 1:length(fnames)) generateReport(fnames[i])
}
    if (any(is.na(fnames)))
        fnames <- choose.files(caption = "Choose Sequence File(s)",
            multi = TRUE, filters = Filters[c("txt", "All"),
                ])
  ampliconReportT(fnames, minMass, maxMass)
  ampliconReportC(fnames, minMass, maxMass)
}