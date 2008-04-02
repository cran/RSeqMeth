`analyzeSequenom` <-
function (matchedName = NA, methName = NA, seqNames = NA, cull=TRUE, 
  writeOut = TRUE, quiet = FALSE, minMass = 1500, maxMass = 7000) 
{
    SNRSum <- function(peaks) {
        Sum <- 0
        for (k in 1:length(peaks)) Sum = Sum + peaks[[k]]$SNR
        return(Sum)
    }
    newPeaks <- function(pTable) {
        newTable = list()
        for (i in 1:nrow(pTable)) if (pTable$Ref.type[i] == "NEW " || 
            pTable$Ref.type[i] == "UNKN") 
            newTable = c(newTable, list(pTable[i, ]))
        return(newTable)
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
    makeCleaveT <- function(seq, pos) {
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
    parseTable <- function(filename) {
        if (!quiet) cat("Loading matched peaks table - may take a few minutes\n")
        table <- read.delim(filename, header = TRUE, colClasses = c("character", 
            "numeric", "numeric", "numeric", "numeric", "numeric", 
            "numeric", "numeric", "character", "character"), 
            skip = 1)
        ind <- indexTable(table)
        p <- (length(ind)/2)
        for (k in 1:p) {
            if (k == 1) {
                temp <- list(Name = table[(ind[k] - 1), 1], Samples = 
                parseAmplicon(table[ind[k]:ind[(k+1)], ]))
                ampTable <- list(temp)
            }
            else {
                temp <- list(Name = table[(ind[(((k - 1) * 2) + 
                  1)] - 1), 1], Samples = (parseAmplicon(table[ind[(((k - 
                  1) * 2) + 1)]:ind[(((k - 1) * 2) + 2)], ])))
                ampTable <- c(ampTable, list(temp))
            }
        }
        if (!quiet) cat("Loaded\n")
        return(ampTable)
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
                  cleaves <- list(makeCleaveT(substr(revcom, 1, 
                    i), 1))
                }
                else cleaves <- c(cleaves, list(makeCleaveT(substr(revcom, 
                  (upto + 1), i), (upto + 1))))
                upto = i
            }
        }
        if (upto != (nchar(revcom))) 
            cleaves <- c(cleaves, list(makeCleaveT(substr(revcom, 
                (upto + 1), nchar(revcom)), (upto + 1))))
        return(cleaves)
    }
    parseSequenceC <- function(filename) {
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
                  cleaves <- list(makeCleaveC(substr(revcom, 1, 
                    i), 1))
                }
                else cleaves <- c(cleaves, list(makeCleaveC(substr(revcom, 
                  (upto + 1), i), (upto + 1))))
                upto = i
            }
        }
        if (upto != (nchar(revcom))) 
            cleaves <- c(cleaves, list(makeCleaveC(substr(revcom, 
                (upto + 1), nchar(revcom)), (upto + 1))))
        return(cleaves)
    }
    parseAmplicon <- function(x) {
        samples <- 0
        skipto <- 0
        for (i in 1:nrow(x)) if (i > skipto - 1) {
            if (length(grep("Plate", x[i, 1])) > 0) {
                flag <- 0
                for (j in (i + 1):nrow(x)) if (length(grep("Plate", 
                  x[j, 1])) > 0) {
                  skipto <- j
                  flag <- 1
                  break
                }
                if (flag == 0) 
                  (skipto <- nrow(x))
                temp <- list(Name = x[i, 1], PeakTable = x[(i + 
                  1):(skipto - 1), ])
                temp = c(temp, New = list(newPeaks(temp$PeakTable)))
                if (samples == 0) 
                  samp <- list(temp)
                else {
                  samp <- c(samp, list(temp))
                }
                samples <- samples + 1
            }
        }
        return(samp)
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
    indexTable <- function(x) {
        index <- vector()
        if (!quiet) cat("---------------------\n")
        upto <- 0
        for (i in 1:nrow(x)) {
            if (substr(x[i, 1], 3, 8) == "Methyl") 
                index = c(index, i)
            if (i > upto) {
                if (!quiet) cat(".")
                upto = upto + (nrow(x)/20)
            }
        }
        if (!quiet) cat("\n")
        return(index)
    }
    loadMethFile <- function(name) {
        temp <- read.delim(file = name, header = F, row.names = 1,stringsAsFactors=F)
        tnames <- temp[1,]
        names(temp) <- tnames
        temp <- temp[2:nrow(temp),]
        AmpNames <- character()
        AmpStarts <- numeric()
        for (i in 1:length(tnames)) {
            slist <- vector()
            for (j in nchar(tnames[i]):1) if (substr(tnames[i], 
                j, j) == "_") 
                slist <- c(slist, j)
            pp <- slist[2]
            if (length(AmpNames) == 0) {
                AmpNames = substr(tnames[i], 1, (pp - 1))
                AmpStarts = i
            }
            else if (!(substr(tnames[i], 1, (pp - 1)) == AmpNames[length(AmpNames)])) {
                AmpNames = c(AmpNames, substr(tnames[i], 1, (pp - 
                  1)))
                AmpStarts = c(AmpStarts, i)
            }
        }
        AmpStarts = c(AmpStarts, (length(tnames) + 1))
        methTable <- list()
        for (i in 1:length(AmpNames)) {
            tempTable <- list(Name = AmpNames[i], Ratios = temp[, 
                AmpStarts[i]:(AmpStarts[i + 1] - 1)])
            if (length(methTable) == 0) 
                methTable <- list(tempTable)
            else methTable <- c(methTable, list(tempTable))
        }
        return(methTable)
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
                  units[[dup]]$FromFrag = i
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites, 
                    "_OL_", temp$CpGSites, sep = "")
                }
                else if (length(cleaves[[i]]$Mass) < length(units[[dup]]$Mass)) {
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites, 
                    "_OL_", cleaves[[i]]$CpGSites, sep = "")
                }
                else {
                  units[[dup]]$CpGSites = paste(units[[dup]]$CpGSites, 
                    "_D_", cleaves[[i]]$CpGSites, sep = "")
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
    calcProportions <- function(samples, CpGs, Cleaves) {
        noUnits <- length(CpGs$Frags)
        colnames <- vector(length = 0, mode = "character")
        for (i in 1:noUnits) {
            for (j in 0:CpGs$Frags[[i]]$CpGs) {
                thename <- paste(CpGs$Frags[[i]]$Name, " ", j, 
                  "/", CpGs$Frags[[i]]$CpGs, sep = "")
                colnames = c(colnames, thename)
            }
        }
        noSamples <- length(samples)
        samplenames <- vector(noSamples, mode = "character")
        for (i in 1:noSamples) samplenames[i] = samples[[i]]$Name
        for (i in 1:noSamples) {
            temp <- calcMethProportion(samples[[i]]$PeakTable, 
                CpGs, Cleaves)
            if (i == 1) 
                result = data.frame(temp)
            else result = data.frame(result, temp)
        }
        names(result) <- samplenames
        row.names(result) <- colnames
        result = data.frame(t(result))
        return(result)
    }
    findFrag <- function(pMass, peaks, tol = 1) {
        for (i in 1:length(peaks)) if ((peaks[[i]]$Mass + tol) > 
            pMass && (peaks[[i]]$Mass - tol) < pMass) 
            return(peaks[[i]])
        return(list(Frag = NA))
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
    calcTableSize <- function(CpGs) {
        ncols <- 0
        for (i in 1:length(CpGs$Frags)) ncols <- ncols + CpGs$Frags[[i]]$CpGs
        return(ncols)
    }
    checkMultipleSilent <- function(fragRefsC, fragRefsA) {
        RefD <- fragRefsA - fragRefsC
        if (length(which(RefD > 0)) > 1) 
            return(FALSE)
        else return(TRUE)
    }
    checkNAs <- function(peakList) {
        naList = vector(length = length(peakList))
        for (k in 1:length(peakList)) naList[k] = peakList[[k]]$SNR
        if (all(is.na(naList))) 
            return(FALSE)
        else return(TRUE)
    }
    checkRefs <- function(fragRefsC) {
        ref <- fragRefsC[1]
        for (k in 2:length(fragRefsC)) if (ref != fragRefsC[k]) 
            return(FALSE)
        return(TRUE)
    }
    cullMethData <- function(mTable, data, unique) {
        mTable <- as.matrix(mTable)
        data <- as.matrix(data)
        mdata <- mat.or.vec(nrow(data), ncol(mTable))
        names(mdata) <- names(mTable)
        rnames <- row.names(data)
        rmnames <- row.names(mTable)
        row.names(mdata) <- rnames
        for (i in 1:nrow(mdata)) {
            for (j in 1:nrow(mTable)) {
                if (substr(rnames[i], 1, (nchar(rnames[i]) - 
                  13)) == rmnames[j]) 
                  break
            }
            mdata[i, ] = mTable[j, ]
        }
        for (i in 1:ncol(data)) for (j in 1:nrow(data)) if (is.na(mdata[j, 
            unique$Frags[[i]]$FromFrag])) 
            data[j, i] = NA
        return(data)
    }
    cullProportionData <- function(mTable, data, unique) {
    	mTable <- as.matrix(mTable)
    	data <- as.matrix(data)
    	mdata <- mat.or.vec(nrow(data), ncol(mTable))
    	names(mdata) <- names(mTable)
    	rnames <- row.names(data)
    	rmnames <- row.names(mTable)
    	row.names(mdata) <- rnames
    	for (i in 1:nrow(mdata)) {
    		for (j in 1:nrow(mTable)) {
    			if (substr(rnames[i],1,(nchar(rnames[i])-13))==rmnames[j]) break
    		}
    		mdata[i,] = mTable[j,]
    	}
    	colnum <- 1
    	for (i in 1:length(unique$Frags)) {
    		noPeaks <- length(unique$Frags[[i]]$Mass)
    		for (j in 1:nrow(data)) if (is.na(mdata[j,unique$Frags[[i]]$FromFrag])) data[j,colnum:(colnum+noPeaks-1)] = NA
    		colnum = colnum + noPeaks
    	}
    	return(data)
    }
    calcSamples <- function(samples, CpGs, Cleaves, weighted = FALSE) {
        noUnits <- length(CpGs$Frags)
        colnames <- vector(length = noUnits, mode = "character")
        for (i in 1:noUnits) colnames[i] = CpGs$Frags[[i]]$Name
        noSamples <- length(samples)
        samplenames <- vector(noSamples, mode = "character")
        for (i in 1:noSamples) samplenames[i] = samples[[i]]$Name
        for (i in 1:noSamples) {
            temp <- calcMeth(samples[[i]]$PeakTable, CpGs, Cleaves, 
                weighted)
            if (i == 1) 
                result = data.frame(temp)
            else result = data.frame(result, temp)
        }
        names(result) <- samplenames
        row.names(result) <- colnames
        result = data.frame(t(result))
        return(result)
    }
    calcPerc <- function(heights) {
        CpGs <- length(heights) - 1
        weights <- (0:CpGs)/CpGs
        wHeights <- heights * weights
        return(round(sum(wHeights)/sum(heights), 2))
    }
    calcMeth <- function(pTable, CpGs, Cleaves, weighted = FALSE) {
        meth <- vector(length = length(CpGs$Frags))
        for (i in 1:length(CpGs$Frags)) {
            noPeaks <- length(CpGs$Frags[[i]]$Mass)
            peaks <- list(length = noPeaks)
            Cfrags <- vector(length = noPeaks)
            Afrags <- vector(length = noPeaks)
            scaledH <- vector(length = noPeaks)
            for (j in 1:noPeaks) {
                peaks[[j]] = findPeak(CpGs$Frags[[i]]$Mass[j], 
                  pTable)
                Cfrags[j] = length(findFrag(CpGs$Frags[[i]]$Mass[j], 
                  CpGs$Peaks)$Frag)
                Afrags[j] = length(findFrag(CpGs$Frags[[i]]$Mass[j], 
                  Cleaves$Peaks)$Frag)
            }
            if (checkNAs(peaks) && checkMultipleSilent(Cfrags, 
                Afrags)) {
                for (j in 1:noPeaks) if (is.na(peaks[[j]]$SNR)) 
                  peaks[[j]]$SNR = 0
                for (j in 1:noPeaks) if ((Afrags[j] - Cfrags[j]) > 
                  0) {
                  scaledH[j] = peaks[[j]]$SNR - ((SNRSum(peaks)/Afrags[j]) * 
                    (Afrags[j] - Cfrags[j]))
                }
                else scaledH[j] = peaks[[j]]$SNR
                if (weighted) 
                  meth[i] = calcPerc(scaledH)
                else meth[i] = round(sum(scaledH[2:noPeaks])/(sum(scaledH)), 
                  2)
                if (is.nan(meth[i])) 
                  meth[i] = NA
                if (!is.na(meth[i])) 
                  if (meth[i] > 1) 
                    meth[i] = 1
            }
            else meth[i:(i + noPeaks)] = NA
        }
        return(meth)
    }
    calcMethProportion <- function(pTable, CpGs, Cleaves) {
        meth <- vector(length = calcTableSize(CpGs))
        column <- 1
        for (i in 1:length(CpGs$Frags)) {
            noPeaks <- length(CpGs$Frags[[i]]$Mass)
            peaks <- list(length = noPeaks)
            Cfrags <- vector(length = noPeaks)
            Afrags <- vector(length = noPeaks)
            scaledH <- vector(length = noPeaks)
            for (j in 1:noPeaks) {
                peaks[[j]] = findPeak(CpGs$Frags[[i]]$Mass[j], 
                  pTable)
                Cfrags[j] = length(findFrag(CpGs$Frags[[i]]$Mass[j], 
                  CpGs$Peaks)$Frag)
                Afrags[j] = length(findFrag(CpGs$Frags[[i]]$Mass[j], 
                  Cleaves$Peaks)$Frag)
            }
            if (checkNAs(peaks) && checkMultipleSilent(Cfrags, 
                Afrags)) {
                for (j in 1:noPeaks) if (is.na(peaks[[j]]$SNR)) 
                  peaks[[j]]$SNR = 0
                for (j in 1:noPeaks) if ((Afrags[j] - Cfrags[j]) > 
                  0) {
                  scaledH[j] = peaks[[j]]$SNR - ((SNRSum(peaks)/Afrags[j]) * 
                    (Afrags[j] - Cfrags[j]))
                }
                else scaledH[j] = peaks[[j]]$SNR
                meth[column:(column + noPeaks - 1)] = scaledH/sum(scaledH)
            }
            else meth[column:(column + noPeaks - 1)] = NA
            column = column + noPeaks
        }
        return(meth)
    }
    
    if (is.na(matchedName))
    matchedName <- choose.files(caption = "Select Matched Peaks Table", 
        multi = FALSE, filters = Filters[c("txt", "All"), ])
    tab <- parseTable(matchedName)
    if (cull) {
          if (is.na(methName))
          methName <- choose.files(caption = "Select Sequenom Methylation Grid", 
            multi = FALSE, filters = Filters[c("txt", "All"), ])
          methTable <- loadMethFile(methName)
    }
    if (any(is.na(seqNames))) {
      seqNames = choose.files(caption = "Choose sequence files to be used in this analysis", 
          multi = TRUE, filters = Filters[c("txt", "All"), ])
      interact <- TRUE
    } else {
      if (length(tab)==length(seqNames)) interact <- FALSE else {
        stop(paste(length(tab),"amplicons in",matchedName,length(seqNames),
          "amplicons in supplied amplicon list"))
      }

    }
    for (i in 1:length(tab)) {
        if (interact)
        AMPa <- parseSequence(select.list(seqNames, title = paste("Sequence file for", 
            tab[[i]]$Name)))
        else AMPa <- parseSequence(seqNames[i])
        AMPc <- findCpGs(AMPa)
        AMPu <- findUniqueUnits(AMPc, minMass, maxMass, skipOL = TRUE)
        AMPuOL <- findUniqueUnits(AMPc, minMass, maxMass, skipOL = FALSE)
        allUnits <- list(Frags = AMPa, Peaks = makePeaks(AMPa))
        CpGUnits <- list(Frags = AMPc, Peaks = makePeaks(AMPc))
        UniqueUnits <- list(Frags = AMPu, Peaks = makePeaks(AMPc))
        UniqueUnitsOL <- list(Frags = AMPuOL, Peaks = makePeaks(AMPc))
        tempUSEQ <- list(Name = tab[[i]]$Name, Ratios = calcSamples(tab[[i]]$Samples, 
            UniqueUnits, allUnits, weighted = FALSE))
        tempUW <- list(Name = tab[[i]]$Name, Ratios = calcSamples(tab[[i]]$Samples, 
            UniqueUnits, allUnits, weighted = TRUE))
        tempUP <- list(Name = tab[[i]]$Name, Ratios = calcProportions(tab[[i]]$Samples, 
            UniqueUnitsOL, allUnits))
        if (cull) {
              tempUSEQ$Ratios = cullMethData(methTable[[i]]$Ratios, 
                  tempUSEQ$Ratios, UniqueUnits)
              tempUW$Ratios = cullMethData(methTable[[i]]$Ratios, tempUW$Ratios, 
                  UniqueUnits)
              tempUP$Ratios <- cullProportionData(methTable[[i]]$Ratios, 
                  tempUP$Ratios, UniqueUnitsOL)
          }
        if (i == 1) {
            dataUSEQ <- list(tempUSEQ)
            dataUW <- list(tempUW)
            dataUP <- list(tempUP)
        }
        else {
            dataUSEQ <- c(dataUSEQ, list(tempUSEQ))
            dataUW <- c(dataUW, list(tempUW))
            dataUP <- c(dataUP, list(tempUP))
        }
    }
    if (writeOut) {
        for (i in nchar(matchedName):1) if (substr(matchedName, i, i) == ".") 
            break
        fname2 <- paste(substr(matchedName, 1, (i - 1)), " Weighted", substr(matchedName, 
            i, nchar(matchedName)), sep = "")
        fname3 <- paste(substr(matchedName, 1, (i - 1)), " Proportion", 
            substr(matchedName, i, nchar(matchedName)), sep = "")
        options(warn = -1)
        for (i in 1:length(dataUW)) {
            if (i == 1) 
                write(file = fname2, dataUW[[i]]$Name, append = FALSE, 
                    sep = "\t")
            else write(file = fname2, dataUW[[i]]$Name, append = TRUE, 
                sep = "\t")
            write.table(file = fname2, dataUW[[i]]$Ratios, append = TRUE, 
                sep = "\t", col.names = NA, row.names = TRUE)
        }
        if (!quiet) cat("Weighted methylation data written to", fname2, "\n")
        for (i in 1:length(dataUP)) {
            if (i == 1) 
                write(file = fname3, dataUP[[i]]$Name, append = FALSE, 
                    sep = "\t")
            else write(file = fname3, dataUP[[i]]$Name, append = TRUE, 
                sep = "\t")
            write.table(file = fname3, dataUP[[i]]$Ratios, append = TRUE, 
                sep = "\t", col.names = NA, row.names = TRUE)
        }
        if (!quiet) cat("Proportion methylation data written to", fname3, "\n")
        options(warn = 0)
    }
    allData <- list("Sequenom" = dataUSEQ, "Weighted" = dataUW, "Proportion" = dataUP)
}


    