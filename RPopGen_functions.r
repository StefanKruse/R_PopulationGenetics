### 		<<<-------------------------------------------->>> ###
### 		<<<---			Taymyr-PopGen Analyses		--->>> ###
###			<<<---			      26.01.2018			--->>> ###
### 		<<<---			Stefan.Kruse@awi.de			--->>> ###
### 		<<<-------------------------------------------->>> ###

# as used in:
# High gene flow and complex treeline dynamics of Larix Mill. stands on the Taymyr Peninsula (north-central Siberia) revealed by nuclear microsatellites
# ... by Kruse, S., Epp, L. S., Wieczorek, M., Pestryakova, L. A., Stoof-Leichsenring, K. R., Herzschuh, U.
# ... in Tree Genetics & Genomes
# ... doi:10.1007/s11295-018-1235-3


print(" --- start loading functions file --- ")

### 		<<<-------------------------------------------->>> ###
### 		<<<---		reformat genind-objects			--->>> ###

			# ... internal function to converse allele-tab to merged allele sizes per locus
			getallels <- function(genindallels, allelmissing="000000")
			{
				if(is.na(sum(genindallels))==TRUE)
				{
					return(allelmissing)
				} else
				{
					if(sum(genindallels)==0)
					{
						return(allelmissing)
					} else
					{
						allelselect=substr(names(genindallels),6,8)[which(genindallels!=0)]
						if(length(allelselect)==1)
						{
							return(paste0(allelselect,allelselect))
						} else
						{
							return(paste0(allelselect[1],allelselect[2]))
						}
					}
				}
			}
			print(" --> getallels(genindallels-object, allelmissing-identifier)")
			genind2table <- function(geninddf)
			{

				# adapt names and save data.names
				namesout=row.names(geninddf$tab)
				namesadapted=paste0(as.character(pop(setPop(geninddf,~Site))),"_",as.character(pop(setPop(geninddf,~AgeClass))))
				namesadapted=make.unique(namesadapted,sep="_")
					
				# convert allel counts to allel colums 
				
						# convert alleles to 6-characters, e.g. "000000"
						msatdfout=NULL
						for(msat in levels(geninddf@loc.fac))
						{

							newallels=apply(geninddf$tab[,geninddf@loc.fac==msat],1,getallels)
							
							# split alleles to 2 cols
							# exchange "000" by "-9"
							newallels=gsub("000000", " -9 -9",newallels)
							
							if(msat==levels(geninddf@loc.fac)[1])
							{
								msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
							} else
							{
								msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
							}
						}
						lvlout=NULL;for(i in 1:length(levels(geninddf@loc.fac)))lvlout=c(lvlout,rep(levels(geninddf@loc.fac)[i],2))
						names(msatdfout)=lvlout
						
						# add names
						msatdfout$Name=namesout
						msatdfout$Name=namesadapted
						msatdfout=cbind(msatdfout, strata(geninddf), other(geninddf))
						
						## sort individuals by populations
							allpopi=as.character(pop(setPop(geninddf,~Site)))
							
							orderi=NULL
							for(popi in allpopsordered[which(allpopsordered%in%unique(allpopi))])
							{
								orderi=c(orderi,which(allpopi==popi))
							}
							msatdfout=msatdfout[orderi,]
						
						# write output
						return(msatdfout)
			}
			print(" --> genind2table(genind-object)")
			genind2tableROH <- function(geninddf)
			{
				# adapt names and save data.names
				namesout=row.names(geninddf$tab)
				popsout=as.character(pop(setPop(geninddf,~Site)))
					
				# convert allel counts to allel colums 
				# convert alleles to 6-characters, e.g. "000000"
						msatdfout=NULL
						for(msat in levels(geninddf@loc.fac))
						{

							newallels=apply(geninddf$tab[,geninddf@loc.fac==msat],1,getallels)
							
							# split alleles to 2 cols
							# exchange "000" by "-9"
							newallels=gsub("000000", " -9 -9",newallels)
							
							if(msat==levels(geninddf@loc.fac)[1])
							{
								msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
							} else
							{
								msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
							}
						}
						lvlout=NULL;for(i in 1:length(levels(geninddf@loc.fac)))lvlout=c(lvlout,rep(levels(geninddf@loc.fac)[i],2))
						names(msatdfout)=lvlout
						
						# add names
						msatdfout=cbind(Name=namesout, Population=popsout, msatdfout, geninddf$strata, other(geninddf))
						
						# write output
						return(msatdfout)

			}
			print(" --> genind2tableROH(genind-object)")
			genind2tableCSV_id_H <- function(geninddf)
			{
				# adapt names and save data.names
				namesout=row.names(geninddf$tab)
				popsout=as.character(pop(setPop(geninddf,~Site)))
					
				# convert allel counts to allel colums 
				# convert alleles to 6-characters, e.g. "000000"
						msatdfout=NULL
						for(msat in levels(geninddf@loc.fac))
						{

							newallels=apply(geninddf$tab[,geninddf@loc.fac==msat],1,getallels)
							
							# split alleles to 2 cols
							# exchange "000" by "-9"
							newallels=gsub("000000", " -9 -9",newallels)
							
							if(msat==levels(geninddf@loc.fac)[1])
							{
								msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
							} else
							{
								msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
							}
						}
						lvlout=NULL;for(i in 1:length(levels(geninddf@loc.fac)))lvlout=c(lvlout,rep(levels(geninddf@loc.fac)[i],2))
						names(msatdfout)=lvlout
						
						# add names
						msatdfout=cbind(Name=namesout, Population=popsout, msatdfout, geninddf$strata, other(geninddf))
						
						# write output
						return(msatdfout)

			}
			print(" --> genind2tableCSV_id_H(genind-object)")

### 		<<<-------------------------------------------->>> ###
### 		<<<---		 general functions 				--->>> ###

### probability matrix comparison
			compareQmatrix <- function(Qmat1, Qmat2)
			{
				# define permutations
				colpermut=permutations(dim(Qmat1)[2],dim(Qmat1)[2],1:dim(Qmat1)[2])

				# calculate differences for each permutation
				dif=NULL
				for(i in 1:dim(colpermut)[2])
				{
					dif=c(dif,sum(abs(Qmat1-Qmat2[,colpermut[i,]])))
				}
				
				return(data.frame(Permut=apply(colpermut,1,function(x)paste0(x,collapse="-")), Diff=dif/dim(colpermut)[2]))
				
			}
			print(" --> compareQmatrix(Qmat1, Qmat2)")
### linkage disequillibrium
			ld2own <- function(geninddf,nloci)
			{
				ival=jval=T2val=dfval=Pvalval=NULL
				for(i in 1:nloci)
				{
					for(j in 1:nloci)
					{
						if(j>i)
						{
							ld2ij=LD2(as.loci(geninddf),c(i,j))
							ival=c(ival, i)
							jval=c(jval, j)
							T2val=c(T2val, ld2ij$T2[1])
							dfval=c(dfval, ld2ij$T2[2])
							Pvalval=c(Pvalval, ld2ij$T2[3])
						} else
						{
							ival=c(ival, i)
							jval=c(jval, j)
							T2val=c(T2val, NA)
							dfval=c(dfval, NA)
							Pvalval=c(Pvalval, NA)
						}
					}
				}
				
				return(data.frame(ival,jval, T2val, dfval, Pvalval))
			}### FUNCTION END
			print(" --> ld2own(genind-object,nloci)")
			# ... missing genotypes are removed prior analysis, NAtoexclude="pairs" => only of compared Loci | "all" => all individuals prior to analysis
			LDpair <- function(genindx,NAtoexclude="pairs")
			{
				nlocgenindx=nLoc(genindx)
				if(NAtoexclude=="all")
				{
					locigenindx=as.loci(missingno(genindx, "geno", quiet=TRUE))
				}
				ival=jval=T2val=dfval=Pvalval=NULL
				for(i in 1:nlocgenindx)
				{
					for(j in 1:nlocgenindx)
					{
						if(j<i)
						{
							if (NAtoexclude=="pairs")
							{
								locigenindx=as.loci(missingno(genindx[loc=c(i,j)], "geno", quiet=TRUE))
								ld2ij=LD2(locigenindx)
							} else if (NAtoexclude=="all")
							{
								ld2ij=LD2(locigenindx,c(i,j))
							}
							ival=c(ival, i)
							jval=c(jval, j)
							T2val=c(T2val, ld2ij$T2[1])
							dfval=c(dfval, ld2ij$T2[2])
							Pvalval=c(Pvalval, ld2ij$T2[3])
						} else
						{
							ival=c(ival, i)
							jval=c(jval, j)
							T2val=c(T2val, NA)
							dfval=c(dfval, NA)
							Pvalval=c(Pvalval, NA)
						}
					}
				}
				return(list(Pval=matrix(Pvalval, nrow=nlocgenindx, ncol=nlocgenindx, byrow=T, dimnames=list(as.vector(levels(genindx@loc.fac)), as.vector(levels(genindx@loc.fac)))), T2val=matrix(T2val, nrow=nlocgenindx, ncol=nlocgenindx, byrow=T, dimnames=list(as.vector(levels(genindx@loc.fac)), as.vector(levels(genindx@loc.fac)))), DFval=matrix(dfval, nrow=nlocgenindx, ncol=nlocgenindx, byrow=T, dimnames=list(as.vector(levels(genindx@loc.fac)), as.vector(levels(genindx@loc.fac)))))
				)
				
			}
			print(" --> LDpair(genind-object,NAtoexclude=\"pairs\")")
			# ... bootstrapped loci random allels (rows shuffled)
			LDpairboot <- function(genindx,NAtoexclude="pairs", bootn=100)
			{
				nlocgenindx=nLoc(genindx)
				
				t2vallist=dfvallist=pvallist=list()
				
				# set up empty matrix
				for(i in 1:bootn)
				{
					t2valmat=matrix(NA, nrow=nlocgenindx, ncol=nlocgenindx, byrow=T, dimnames=list(as.vector(levels(genindx@loc.fac)), as.vector(levels(genindx@loc.fac))))
					t2vallist[[i]]=t2valmat
					
					dfvalmat=matrix(NA, nrow=nlocgenindx, ncol=nlocgenindx, byrow=T, dimnames=list(as.vector(levels(genindx@loc.fac)), as.vector(levels(genindx@loc.fac))))
					dfvallist[[i]]=dfvalmat
					
					pvalmat=matrix(NA, nrow=nlocgenindx, ncol=nlocgenindx, byrow=T, dimnames=list(as.vector(levels(genindx@loc.fac)), as.vector(levels(genindx@loc.fac))))
					pvallist[[i]]=pvalmat
				}
				
				

				for(i in 1:nlocgenindx)
				{
					for(j in 1:nlocgenindx)
					{

						if(j>i)
						{
							if (NAtoexclude=="pairs")
							{
								# print progress to console
									print(paste("i=",i,"  -  ","j=",j))
									

								locigenindx=as.loci(missingno(genindx[loc=c(j,i)], "geno", quiet=TRUE))
								for(booti in 1:bootn)
								{

									# [,1]==pop
									# [,2]==Locus1
									# [,3]==Locus2
									locigenindx[,2]=locigenindx[sample(1:dim(locigenindx)[1]),2]
									ld2ij=LD2(locigenindx)
									t2vallist[[booti]][j,i]=ld2ij$T2[1]
									dfvallist[[booti]][j,i]=ld2ij$T2[2]
									pvallist[[booti]][j,i]=ld2ij$T2[3]
								}
								
						
								
							} 
						}
					}
				}
				
				
				return(	list(Pval=pvallist,T2val=t2vallist,DFval=dfvallist) )
			}
			print(" --> LDpairboot(genind-object,NAtoexclude=\"pairs\",bootn)")

### plot alleles per loci
			#... allele per locus plot per strata in colors
			plotLociPerMSat <- function(genindx, freq=TRUE, plotting=TRUE)
			{
				nloci=nLoc(genindx)
				sqposi=ceiling(sqrt(nloci))
				popi=levels(pop(genindx))
				cols=rainbow(length(popi))
				
				if(plotting)dev.new(height=1.5*sqposi, width=1.5*sqposi)
				if(plotting)par(mfrow=c(sqposi,sqposi), mar=c(2,2,2,0))
					allelsout=list()
				for(loclvli in levels(genindx$loc.fac))
				{
					allelsii=NULL
					for(popii in popi)
					{
						valsi=apply(popsub(genindx,popii)$tab[,which(popsub(genindx,popii)$loc.fac==loclvli)],2,function(x)sum(na.omit(x)))
						if(freq)valsi=valsi/sum(valsi)
						allelsi=as.numeric(unlist(strsplit(names(valsi),"\\."))[c(FALSE,TRUE)])
						if(which(popi==popii)==1)
						{
							if(freq)
							{
								ylimimax=1
							} else
							{
								ylimimax=max(valsi)*1.1
							}

							if(plotting)plot(valsi~c(allelsi+0.4*(which(popi==popii)-1)),typ="h",lwd=3, main=loclvli, col=cols[which(popi==popii)], ylim=c(0,ylimimax))
						} else
						{
							if(plotting)points(valsi~c(allelsi+0.4*(which(popi==popii)-1)),typ="h",lwd=3, main=loclvli, col=cols[which(popi==popii)])
						}
						
						allelsii=c(allelsii, allelsi)
					}
					allelsout[[which(levels(genindx$loc.fac)==loclvli)]]=allelsii
				}
				legend("topleft", popi, lty=1, lwd=2, col=cols)
				
				names(allelsout)=levels(genindx$loc.fac)
				
				return(allelsout)
			}
			print(" --> plotLociPerMSat(genind-object, freq, plotting)")

### extract alleles per loci
			extractNallelesPerMSat <- function(genindx)
			{
				nallels=NULL
				for(loclvli in levels(genindx$loc.fac))
				{
					valsi=apply(genindx$tab[,which(genindx$loc.fac==loclvli)],2,function(x)sum(na.omit(x)))
					nallels=c(nallels,length(valsi[valsi!=0]))
				}
				return(nallels)
				
			}
			print(" --> extractNallelesPerMSat(genind-object)")

### calculate the number of rare alleles per population
			# ... for genind object
			rareallelscalc_genind=function(genindi)
			{
				valsi=Pop=NULL
				allallfreq=tab(genindi, freq=TRUE, NA.method="zero")
				for(popi in levels(pop(genindi)))
				{
					pospopi=which(pop(genindi)==popi)
					valsi=rbind(valsi,cbind(as.data.frame(t(apply(allallfreq[pospopi,],2,sum))), Pop=popi))
				}
				# calculate frequencies for each population and locus
					valsi=rbind(valsi, cbind(as.data.frame(t(apply(allallfreq,2,sum))), Pop="TOTAL"))
					valsfreq=cbind(valsi[,-which(names(valsi)=="Pop")]/apply(valsi[,-which(names(valsi)=="Pop")],1,sum), Pop=valsi[,which(names(valsi)=="Pop")])

				# calculate frequencies for each locus
					freqpopout_overall=data.frame(Pop=valsi$Pop)
					for(loclvli in levels(genindi$loc.fac))
					{
						subvalsi=valsi[,which(substr(names(valsi),1,4)==loclvli)]
						freqpopout_overall=cbind(freqpopout_overall,subvalsi/apply(subvalsi,1,sum))
					}

					# .. select  <1%
					smaller1=names(freqpopout_overall[freqpopout_overall$Pop=="TOTAL",-which(names(freqpopout_overall)=="Pop")])[which(freqpopout_overall[freqpopout_overall$Pop=="TOTAL",-which(names(freqpopout_overall)=="Pop")]<0.01)]
					sub=freqpopout_overall[,c("Pop",smaller1)]

				# count number of rare alleles
					Nrareallel=NULL
					for(popi in levels(pop(genindi)))
					{
						Nrarealleli=locnami=NULL
						for(loclvli in levels(genindi$loc.fac))
						{
							if(length(which(substr(names(sub),1,4)==loclvli))>0)
							{
								plocii=which(substr(names(sub),1,4)==loclvli)
								Nrarealleli=c(Nrarealleli,length(which(sub[sub$Pop==popi,plocii]!=0)))
							} else
							{
								Nrarealleli=c(Nrarealleli,NA)
							}
							locnami=c(locnami, loclvli)
						}
						
						Nrareallel=rbind(Nrareallel,data.frame(Loc=locnami, NRA=Nrarealleli, Pop=popi))
					}
					nraout=with(Nrareallel	, aggregate(NRA, list(Pop), sum))
					names(nraout)=c("Pop", "NRA")
					
				return(nraout)
			}
			print(" --> rareallelscalc_genind(genind-object)")
			# ... for list object
			areallelscalc=function(genindi)
			{
				valsi=Pop=NULL
				allallfreq=tab(genindi, freq=TRUE, NA.method="zero")
				for(popi in levels(pop(genindi)))
				{
					pospopi=which(pop(genindi)==popi)
					valsi=rbind(valsi,cbind(as.data.frame(t(apply(allallfreq[pospopi,],2,sum))), Pop=popi))
				}

				# calculate frequences for each population and locus
				valsi=rbind(valsi, cbind(as.data.frame(t(apply(allallfreq,2,sum))), Pop="TOTAL"))

					# ... calculate frequences for each locus 
					freqpopout_overall=data.frame(Pop=valsi$Pop)
					for(loclvli in levels(genindi$loc.fac))
					{
						subvalsi=valsi[,which(substr(names(valsi),1,4)==loclvli)]
						freqpopout_overall=cbind(freqpopout_overall,subvalsi/apply(subvalsi,1,sum))
					}
					
					# .. seelct  <1%
					smaller1=names(freqpopout_overall[freqpopout_overall$Pop=="TOTAL",-which(names(freqpopout_overall)=="Pop")])[which(freqpopout_overall[freqpopout_overall$Pop=="TOTAL",-which(names(freqpopout_overall)=="Pop")]<0.01)]
					sub=freqpopout_overall[,c("Pop",smaller1)]

						# count for each locus populations !=0
						Nrareallel=NULL
						for(popi in levels(pop(genindi)))
						{
							Nrarealleli=locnami=NULL
							for(loclvli in levels(genindi$loc.fac))
							{
								if(length(which(substr(names(sub),1,4)==loclvli))>0)
								{
									plocii=which(substr(names(sub),1,4)==loclvli)
									Nrarealleli=c(Nrarealleli,length(which(sub[sub$Pop==popi,plocii]!=0)))
								} else
								{
									Nrarealleli=c(Nrarealleli,NA)
								}
								locnami=c(locnami, loclvli)
							}
							
							Nrareallel=rbind(Nrareallel,data.frame(Loc=locnami, NRA=Nrarealleli, Pop=popi))
						}
						nraout=with(Nrareallel	, aggregate(NRA, list(Pop), sum))
						names(nraout)=c("Pop", "NRA")
						
						return(nraout)
			}
			print(" --> areallelscalc(genind-object)")

	
### m test for correlations
			cor.mtest <- function(mat, conf.level = 0.95)
			{
				mat <- as.matrix(mat)
				n <- ncol(mat)
				p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
				diag(p.mat) <- 0
				diag(lowCI.mat) <- diag(uppCI.mat) <- 1
				for(i in 1:(n-1)){
					for(j in (i+1):n){
						tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
						p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
						lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
						uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
					}
				}
				return(list(p.mat, lowCI.mat, uppCI.mat))
			}
			print(" --> cor.mtest(matrix, conf.level)")

### plot pairwise significant different group labels
			# ... function to add group labels to not significant different groups
			# ... adapted from https://www.r-bloggers.com/automated-determination-of-distribution-groupings-a-stackoverflow-collaboration/
			plotWilcoxGroups <- function(datf, xvalnam="Site", yvalnam, pniv, padjmeth="bonferroni",plottin=TRUE)
			{
				x=datf[,xvalnam]
				y=datf[,yvalnam]
				m1=cbind.data.frame(x,y)

				mw <- pairwise.wilcox.test(m1[,2], m1[,1], p.adj=padjmeth)
				mw$p.value

				# matrix showing connections between levels
				g=mw$p.value
				g=cbind(rbind(0, g), 0)
				g=replace(g, is.na(g), FALSE)
				g=g + t(g)
				g
				diag(g)=1
				rownames(g)=colnames(g)= levels(datf[,xvalnam])
				g=as.matrix(g > pniv)
				g=replace(g, is.na(g), FALSE)
				g=abs(g)
				
				# re-arrange data into an "edge list" for use in igraph (i.e. which groups are "connected") - Solution from "David Eisenstat"
				same=which(g==1)
				k=arrayInd(same, dim(g))

				g2=data.frame(rownames(g)[k[,1]], colnames(g)[k[,2]])
				g2=g2[order(g2[[1]]),]
				g3=simplify(graph.data.frame(g2,directed = FALSE))
				get.data.frame(g3)

				V(g3)$color="white"
				V(g3)$label.color="black"
				V(g3)$size=20

				n=length(levels(datf[,xvalnam]))
				g4=data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
				g4<- g4[order(g4[[1]]),]
				g5=simplify(graph.data.frame(g4,directed = FALSE))

				# calcuate the maximal cliques - these are groupings where every node is connected to all others
				cliq=maximal.cliques(g5) # Solution from "majom"
				cliq2=lapply(cliq, as.numeric)
				# Reorder by level order - Solution from "MrFlick"
				ml<-max(sapply(cliq, length))
				reord=do.call(order, data.frame(
				  do.call(rbind,
						  lapply(cliq2, function(x) c(sort(x), rep.int(0, ml-length(x))))
				  )
				))
				cliq=cliq2[reord]


				# generate labels to factor levels
				lab.txt=vector(mode="list", n)
				lab=letters[seq(cliq)]
				for(i in seq(cliq)){
				  for(j in cliq[[i]]){
					lab.txt[[j]]=paste0(lab.txt[[j]], lab[i])
				  }
				}
				
				if(plottin)
				{	
					text(y=1:length(bp$names), x=bp$stats[4,]+0.01*(mean(na.omit(bp$stats[4,]))), labels=lab.txt, col="black", cex=1, font=2, adj=0)
				}
				
				outtdf=data.frame(Group=unlist(lab.txt)	, Site=levels(datf[,xvalnam]))
				return(outtdf[order(outtdf$Group),])
			}
			print(" --> plotWilcoxGroups(dataframe, xvalnam, yvalnam, pniv, padjmeth, plottin)")

### PCAsignificance function
# ... based on BiodiversityR's PCAsiginifcance()
			DAPCsignificance <- function (dapc, axes = 8) 
			{
				eigen <- dapc$eig
				tot <- sum(eigen)
				p <- length(eigen)
				if (p < axes) 
				{
					axes <- p
				}
				varexplained <- array(dim = c(7, p))
				varexplained[1, ] <- eigen[1:p]
				varexplained[2, ] <- varexplained[1, ]/tot * 100
				varexplained[3, 1] <- varexplained[2, 1]
				for (i in 2:p) {
					varexplained[3, i] <- varexplained[3, i - 1] + varexplained[2, 
						i]
				}
				for (i in 1:p) {
					varexplained[6, i] <- 1/i
				}
				for (i in 1:p) {
					varexplained[4, i] <- sum(varexplained[6, i:p])/p * 100
				}
				varexplained[5, 1] <- varexplained[4, 1]
				for (i in 2:p) {
					varexplained[5, i] <- varexplained[5, i - 1] + varexplained[4, 
						i]
				}
				for (i in 1:p) {
					if (varexplained[2, i] > varexplained[4, i]) {
						varexplained[6, i] <- TRUE
					}
					else {
						varexplained[6, i] <- FALSE
					}
					if (varexplained[3, i] > varexplained[5, i]) {
						varexplained[7, i] <- TRUE
					}
					else {
						varexplained[7, i] <- FALSE
					}
				}
				rownames(varexplained) <- c("eigenvalue", "percentage of variance", 
					"cumulative percentage of variance", "broken-stick percentage", 
					"broken-stick cumulative %", "% > bs%", "cum% > bs cum%")
				colnames(varexplained) <- c(1:p)
				return(varexplained[, 1:axes])
			}
			print(" --> DAPCsignificance(dapc,axes-number)")

#### functions for analysing list output
			listsummary <- function(listwithtable, pvalout=FALSE)
			{
				mat=matrix(rep(NA, length(listwithtable)*length(as.vector(listwithtable[[1]]))), ncol=length(listwithtable))
				for(repeati in 1:length(listwithtable))
				{
					# ... convert table to vector
					asvec=as.vector(listwithtable[[repeati]])#1. column, alle zeilen, 2. column ....
					# ... buffer colnames/rownames
					colnms=colnames(listwithtable[[repeati]])
					rownms=row.names(listwithtable[[repeati]])
					# ... create matrix of 1:x elements (in columns repeats)
					mat[,repeati]=asvec
				}
				
				nummerge=paste(round(apply(mat, 1, mean),3), round(apply(mat, 1, sd),3), sep="+/-")
				ttestlist=NULL
				for(rowi in 1:dim(mat)[1])
				{
					if(length(na.omit(mat[rowi,]))==0)
					{
						ttestlist=c(ttestlist,NA)
					} else if(length(unique(mat[rowi,]))==1)
					{
						ttestlist=c(ttestlist,NA)
						print(paste("   - t-test not possible because data constant in row: ",rowi))
					} else
					{
						ttestlist=c(ttestlist,t.test(mat[rowi,])$p.value)
					}
				}
				
				pvalvec=ttestlist
				nummerge_pval=formatC(pvalvec,format="g",digits=2)

				xdf4sum=as.data.frame(matrix(apply(mat, 1, mean), nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)))
					xdf4sumpvalue=as.data.frame(matrix(pvalvec, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)))
						tot11=paste(round(apply(xdf4sum[row.names(xdf4sum)%in%loc_11,],2,mean),3),round(apply(xdf4sum[row.names(xdf4sum)%in%loc_11,],2,sd),3),sep="+/-")
						tot9=paste(round(apply(xdf4sum[row.names(xdf4sum)%in%loc_9,],2,mean),3),round(apply(xdf4sum[row.names(xdf4sum)%in%loc_9,],2,sd),3),sep="+/-")
						tot11_pval=paste(round(apply(xdf4sumpvalue,2,mean),3), round(apply(xdf4sumpvalue,2,sd),3), sep="+/-")
						tot9_pval=paste(round(apply(xdf4sumpvalue[row.names(xdf4sumpvalue)%in%loc_9,],2,mean),3),round(apply(xdf4sumpvalue[row.names(xdf4sumpvalue)%in%loc_9,],2,sd),3),sep="+/-")
					outdf=rbind(as.data.frame(matrix(nummerge, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE), tot11, tot9)
					outdf_pval=rbind(as.data.frame(matrix(nummerge_pval, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE), tot11_pval, tot9_pval)
					row.names(outdf)[(dim(outdf)[1]-1):dim(outdf)[1]]=c("Mean11Loci", "Mean9Loci")
					row.names(outdf_pval)[(dim(outdf_pval)[1]-1):dim(outdf_pval)[1]]=c("Mean11Loci", "Mean9Loci")
					if(pvalout)
					{
						return(list(meansd=outdf,pvalues=outdf_pval))
					} else
					{
						return(outdf)
					}
			}
			print(" --> listsummary(listwithtable,pvalout)")
			# ... for 8 loci
			listsummary_8 <- function(listwithtable, pvalout=FALSE)
			{
				mat=matrix(rep(NA, length(listwithtable)*length(as.vector(listwithtable[[1]]))), ncol=length(listwithtable))
				for(repeati in 1:length(listwithtable))
				{
					# ... convert table to vector
					asvec=as.vector(listwithtable[[repeati]])
					# ... buffer colnames/rownames
					colnms=colnames(listwithtable[[repeati]])
					rownms=row.names(listwithtable[[repeati]])
					# ... create matrix of 1:x elements (in columns repeats)
					mat[,repeati]=asvec
				}
				
				nummerge=paste(round(apply(mat, 1, mean),3), round(apply(mat, 1, sd),3), sep="+/-")
				ttestlist=NULL
				for(rowi in 1:dim(mat)[1])
				{
					if(length(na.omit(mat[rowi,]))==0)
					{
						ttestlist=c(ttestlist,NA)
					} else if(length(unique(mat[rowi,]))==1)
					{
						ttestlist=c(ttestlist,NA)
						print(paste("   - ttest nicht moeglich da daten konstant in Zeile:",rowi))
					} else
					{
						ttestlist=c(ttestlist,t.test(mat[rowi,])$p.value)
					}
				}
				pvalvec=ttestlist
				nummerge_pval=formatC(pvalvec,format="g",digits=2)

				xdf4sum=as.data.frame(matrix(apply(mat, 1, mean), nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)))
					xdf4sumpvalue=as.data.frame(matrix(pvalvec, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)))
						tot8=paste(round(apply(xdf4sum,2,mean),3),round(apply(xdf4sum,2,sd),3),sep="+/-")
						tot8_pval=paste(round(apply(xdf4sumpvalue,2,mean),3), round(apply(xdf4sumpvalue,2,sd),3), sep="+/-")
					outdf=rbind(as.data.frame(matrix(nummerge, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE), tot8)
					outdf_pval=rbind(as.data.frame(matrix(nummerge_pval, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE), tot8_pval)
					row.names(outdf)[dim(outdf)[1]]=c("MeanLoci")
					row.names(outdf_pval)[dim(outdf_pval)[1]]=c("MeanLoci")
					if(pvalout)
					{
						return(list(meansd=outdf,pvalues=outdf_pval))
					} else
					{
						return(outdf)
					}
			}
			print(" --> listsummary_8(listwithtable,pvalout)")
			# ... for data frame output
			listsummary4dfs <- function(listwithtable)
			{
				mat=matrix(rep(NA, length(listwithtable)*length(as.vector(unlist(listwithtable[[1]])))), ncol=length(listwithtable))
				for(repeati in 1:length(listwithtable))
				{
					# ... convert table to vector
					asvec=as.vector(unlist(listwithtable[[repeati]]))
					# ... buffer colnames/rownames
					colnms=colnames(listwithtable[[repeati]])
					rownms=row.names(listwithtable[[repeati]])
					# ... create matrix of 1:x elements (in columns repeats)
					mat[,repeati]=asvec
				}
				
				nummerge=paste(round(apply(mat, 1, mean),3), round(apply(mat, 1, sd),3), sep="+/-")
					
				xdf4sum=as.data.frame(matrix(apply(mat, 1, mean), nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)))
					tot11=paste(round(apply(xdf4sum,2,mean),3), round(apply(xdf4sum,2,sd),3), sep="+/-")
					tot9=paste(round(apply(xdf4sum[row.names(xdf4sum)%in%loc_9,],2,mean),3),round(apply(xdf4sum[row.names(xdf4sum)%in%loc_9,],2,sd),3),sep="+/-")
				outdf=as.data.frame(matrix(nummerge, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE)	
				outdf=rbind(outdf, tot11, tot9)
				row.names(outdf)[(dim(outdf)[1]-1):dim(outdf)[1]]=c("Mean11Loci", "Mean9Loci")
				
				return(outdf)
			}
			print(" --> listsummary4dfs(listwithtable)")
			# ... with booted confidence intervals
			listsummary4dfsCIbootfst <- function(listwithtable)
			{
				mat=matrix(rep(NA, length(listwithtable)*length(as.vector(unlist(listwithtable[[1]])))), ncol=length(listwithtable))
				for(repeati in 1:length(listwithtable))
				{
					# ... convert table to vector
					asvec=as.vector(unlist(listwithtable[[repeati]]))
					# ... buffer colnames/rownames
					colnms=colnames(listwithtable[[repeati]])
					rownms=row.names(listwithtable[[repeati]])
					# ... create matrix of 1:x elements (in columns repeats)
					mat[,repeati]=asvec
				}
				
				nummerge_lowerci=apply(mat, 1, mean)-1.96*apply(mat, 1, sd)
				nummerge_upperci=apply(mat, 1, mean)+1.96*apply(mat, 1, sd)
				
				outdf_lower=as.data.frame(matrix(nummerge_lowerci, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE)	
				outdf_upper=as.data.frame(matrix(nummerge_upperci, nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE)	
					
				return(list(lowerCI=outdf_lower,upperCI=outdf_upper))
			}
			print(" --> listsummary4dfsCIbootfst(listwithtable)")
			# ... for calculating mean values
			listsummarymean <- function(listwithtable)
			{
				mat=matrix(rep(NA, length(listwithtable)*length(as.vector(listwithtable[[1]]))), ncol=length(listwithtable))
				for(repeati in 1:length(listwithtable))
				{
					# ... convert table to vector
					asvec=as.vector(listwithtable[[repeati]])
					# ... buffer colnames/rownames
					colnms=colnames(listwithtable[[repeati]])
					rownms=row.names(listwithtable[[repeati]])
					# ... create matrix of 1:x elements (in columns repeats)
					mat[,repeati]=asvec
				}
				
				return(as.data.frame(matrix(apply(mat, 1, mean), nrow=length(rownms), ncol=length(colnms), dimnames=list(rownms, colnms)), stringsAsFactors=FALSE)	)
			}
			print(" --> listsummarymean(listwithtable)")



			
			
			
# ... program specific export functions 	
# ... ... GENEPOP
			genind2genepop <- function(xgenind, header, popstrata)
			{
				titlerow=header
				locusnamesrow=paste(levels(xgenind@loc.fac), collapse=",")

				# reformat allele data
					xtab=as.data.frame(xgenind@tab)
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
						# print(msat)
						
						newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallels)
						
						if(msat==levels(xgenind@loc.fac)[1])
						{#neu
							msatdfout=data.frame(newallels)
						} else
						{
							msatdfout=cbind(msatdfout, newallels)
						}
					}
					names(msatdfout)=levels(xgenind@loc.fac)

				# for each population
					popname=paste(as.character(pop(xgenind)),gsub("\\.","_",gsub(",","_",gsub(" ", "_", row.names(msatdfout)))), sep=" ")
					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
					popnameallels=paste(popname, allelsi, sep=" , ")
					
					outputlines=rbind(titlerow, locusnamesrow)
					for(popi in levels(pop(xgenind)))
					{
						outputlines=rbind(outputlines, "Pop")
						for(rowi in which(as.character(pop(xgenind))==popi))
						{
							outputlines=rbind(outputlines, popnameallels[rowi])
						}
					}
					
					return(outputlines)
			}
			print(" --> genind2genepop(genind-object, header, popstrata)")
# ... ... BOTTLENECK
			genind2genepop4BOTTLENECK <- function(xgenind, header)
			{
				titlerow=header
				locusnamesrow=NULL
				for(i in 1:length(levels(xgenind@loc.fac)))
				{
					locusnamesrow=rbind(locusnamesrow,levels(xgenind@loc.fac)[i])
				}
				
				# reformat allele data
					xtab=as.data.frame(xgenind@tab)
					
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
						newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallels)
						
						if(msat==levels(xgenind@loc.fac)[1])
						{
							msatdfout=data.frame(newallels)
						} else
						{
							msatdfout=cbind(msatdfout, newallels)
						}
					}
					names(msatdfout)=levels(xgenind@loc.fac)

				# for each population
					popname=paste(as.character(pop(xgenind)),gsub("\\.","_",gsub(",","_",gsub(" ", "_", row.names(msatdfout)))), sep="ind")
					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
					popnameallels=paste(popname, allelsi, sep=" , ")
					
					outputlines=rbind(titlerow, locusnamesrow)
					for(popi in levels(pop(xgenind)))
					{
						outputlines=rbind(outputlines, "Pop")
						for(rowi in which(as.character(pop(xgenind))==popi))
						{
							outputlines=rbind(outputlines, popnameallels[rowi])
						}
					}
					
				return(outputlines)
					
			}
			print(" --> genind2genepop4BOTTLENECK(genind-object, header)")
# ... ... FreeNA
			genind2genepop4FreeNA <- function(xgenind, header)
			{
				titlerow=header
				locusnamesrow=NULL
				for(i in 1:length(levels(xgenind@loc.fac)))
				{
					locusnamesrow=rbind(locusnamesrow,levels(xgenind@loc.fac)[i])
				}
				
				# reformat allele data
					xtab=as.data.frame(xgenind@tab)

					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
						newallels=apply(xtab[,xgenind@loc.fac==msat],1,function(x)getallels(genindallels=x,allelmissing="999999"))
						
						if(msat==levels(xgenind@loc.fac)[1])
						{
							msatdfout=data.frame(newallels)
						} else
						{
							msatdfout=cbind(msatdfout, newallels)
						}
					}
					names(msatdfout)=levels(xgenind@loc.fac)

					# for each population
						popname=paste(as.character(pop(xgenind)),gsub("\\.","_",gsub(",","_",gsub(" ", "_", row.names(msatdfout)))), sep="ind")
						allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
						popnameallels=paste(popname, allelsi, sep=" , ")
						
						outputlines=rbind(titlerow, locusnamesrow)
						for(popi in levels(pop(xgenind)))
						{
							outputlines=rbind(outputlines, "Pop")
							for(rowi in which(as.character(pop(xgenind))==popi))
							{
								outputlines=rbind(outputlines, popnameallels[rowi])
							}
						}
					
					return(outputlines)
					
			}
			print(" --> genind2genepop4FreeNA(genind-object, header)")
# ... ... INEST2							
			genind2genepop4INEST2 <- function(xgenind,allorpops="all",filebase)
			{
				locusnamesrow=NULL
				for(i in 1:length(levels(xgenind@loc.fac)))
				{
					locusnamesrow=rbind(locusnamesrow,levels(xgenind@loc.fac)[i])
				}

				# reformat allele data
					xtab=as.data.frame(xgenind@tab)
					
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
						# print(msat)
						
						newallels=apply(xtab[,xgenind@loc.fac==msat],1,function(x)getallels(genindallels=x,allelmissing="000000"))
						
						if(msat==levels(xgenind@loc.fac)[1])
						{
							msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
						} else
						{
							msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
						}
					}
					lvlout=NULL;for(i in 1:length(levels(xgenind@loc.fac)))lvlout=c(lvlout,rep(levels(xgenind@loc.fac)[i],2))
					names(msatdfout)=lvlout

				# for each population
					indname=1:nInd(xgenind)
					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse="\t")))
					popnameallels=paste(indname, allelsi, sep="\t")
					
					if(allorpops=="all")
					{
						titlerow=paste0(nInd(xgenind),"\t",nLoc(xgenind), "\t","0")
						outputlines=c(titlerow, locusnamesrow, popnameallels)
						
						# write file 
						write.table(outputlines,paste0(filebase,".txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)
						
					} else if(allorpops=="pops")## not functional yet
					{
						for(popi in levels(pop(xgenind)))
						{
							titlerow=paste0(length(which(as.character(pop(xgenind))==popi)),"\t",nLoc(xgenind), "\t","0")
							outputlinespop=c(titlerow, locusnamesrow)
							for(rowi in which(as.character(pop(xgenind))==popi))
							{
								outputlinespop=c(outputlinespop, popnameallels[rowi])
							}
							# write file
							write.table(outputlinespop,paste0(filebase,"_pop_",popi,".txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)

						}
					}
				
			}
			print(" --> genind2genepop4INEST2(genind-object, allorpops,filebase)")
# ... ... Arlequin
			genind2genepop4Arlequin <- function(xgenind, header, popstrata)
			{
				titlerow=header
				locusnamesrow=paste(levels(xgenind@loc.fac), collapse=",")

				# reformat allele data
					xtab=as.data.frame(xgenind@tab)

					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
						newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallels)
						
						if(msat==levels(xgenind@loc.fac)[1])
						{
							msatdfout=data.frame(newallels)
						} else
						{
							msatdfout=cbind(msatdfout, newallels)
						}
					}
					names(msatdfout)=levels(xgenind@loc.fac)

				# for each population
					popname=as.character(pop(xgenind))
					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
					popnameallels=paste(popname, allelsi, sep=" , ")
					
					outputlines=rbind(titlerow, locusnamesrow)
					for(popi in levels(pop(xgenind)))
					{
						outputlines=rbind(outputlines, "Pop")
						for(rowi in which(as.character(pop(xgenind))==popi))
						{
							outputlines=rbind(outputlines, popnameallels[rowi])
						}
					}
					
				return(outputlines)
					
			}
			print(" --> genind2genepop4Arlequin(genind-object, header, popstrata)")
			genind2arlequin <- function(geninddf,outputdir, filename)
			{
				# format conventions
					# 2 rows per indi
					# 1st col= numbers increasing for each individual
					# ... save original names in extra file
					# 2nd col= population origin
					# 1:12 col = allels loci1:12
					# missing = -9

				header=rbind(
				"[Profile]",
				"Title=\"Microsatellite Larix Taimyr\"",
				paste0("NbSamples=",length(levels(pop(geninddf)))),
				"DataType=MICROSAT",
				"GenotypicData=1",
				"GameticPhase=0",
				"MissingData=\"?\"",
				"LocusSeparator=WHITESPACE",
				paste0("#Loci(N=",nLoc(geninddf),"): ",paste(levels(geninddf@loc.fac), collapse=", "))
				, " "
				, "[Data]"
				, " "
				,"[[Samples]]"
				)
				
				# individuals:
					# Name(unique) SPACE 1 SPACE loc1/1 SPACE loc2/1 ...
					#                      SPACE loc1/2 SPACE loc2/2 ...

					# convert allele counts to allele columns 
						# convert alleles to 6-characters, e.g. "000000"
						msatdfout=NULL
						for(msat in levels(geninddf@loc.fac))
						{

							newallels=apply(geninddf@tab[,geninddf@loc.fac==msat],1,getallels)
							
							# split alleles to 2 cols
							# replace "000" by "?"
							newallels=gsub("000000", "?  ?  ",newallels)
							
							if(msat==levels(geninddf@loc.fac)[1])
							{
								msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
							} else
							{
								msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
							}
						}
						
						namesout=NULL
						for(lvli in levels(geninddf@loc.fac))
						{
							namesout=c(namesout,paste0(lvli,c(",1",",2")))
						}
						names(msatdfout)=namesout
						
						# convert 1-row to 2-row data for each individuals
						msatdfout$indinr=1:dim(msatdfout)[1]
							x1=msatdfout[,c(paste0(levels(geninddf@loc.fac),",1"),"indinr")]
							x2=msatdfout[,c(paste0(levels(geninddf@loc.fac),",2"),"indinr")]
							names(x2)=names(x1)
							msatoutmerged=rbind(x1,x2)
							msatoutmerged=msatoutmerged[order(msatoutmerged$indinr),]
				
					# add pop
						msatoutmerged$PopulationID=
						iddf=rbind(data.frame(id=1:dim(msatdfout)[1], PopulationID=pop(geninddf)),data.frame(id=1:dim(msatdfout)[1], PopulationID=pop(geninddf)))
						msatoutmerged$PopulationID=iddf[order(iddf$id),]$PopulationID
						names(msatoutmerged)[which(names(msatoutmerged)=="indinr")]="IndividualID"
						msatoutmerged=msatoutmerged[,c((dim(msatoutmerged)[2]-1):(dim(msatoutmerged)[2]), 1:(dim(msatoutmerged)[2]-2))]
						names(msatoutmerged)=gsub(",1", "", names(msatoutmerged))
						
						
					# for each population
						poprows=NULL
						for(lvli in levels(pop(geninddf)))
						{
							posindi=which(pop(geninddf)==lvli)
							nindi=length(posindi)
							
							poprows=rbind(poprows,
								paste0("SampleName=\"",lvli,"\" "),
								paste0("SampleSize=",nindi),
								"SampleData={ "
							)
							
							for(indii in posindi)
							{
								poprows=rbind(poprows,
									paste(formatC(indii, width=4), "1",paste(as.matrix(msatoutmerged[msatoutmerged$IndividualID%in%indii,levels(geninddf@loc.fac)][1,]), collapse=" "), collapse=" "),
									paste(formatC("", width=4), " ",paste(as.matrix(msatoutmerged[msatoutmerged$IndividualID%in%indii,levels(geninddf@loc.fac)][1,]), collapse=" "), collapse=" ")
									)
							}
							
							poprows=rbind(poprows, "} ")
						}
						
						
					structrows=rbind(
							"[[Structure]]",
							" ",
							"StructureName = \"Test\"",
							paste0("NbGroups=",1),
							" ",
							"Group =  { ",
							as.matrix(levels(pop(geninddf))),
							"} "
						)
							
						
						
					# output
						write.table(rbind(header,poprows,structrows), paste0(outputdir,filename), col.names=FALSE, row.names=FALSE,, quote=FALSE)
							
			}
			print(" --> genind2arlequin(genind-object,outputdir,filename)")
# ... ... diveRsity
			genind2genepop4diveRsity <- function(xgenind, header, popstrata)
			{
				titlerow=header
				locusnamesrow=paste(levels(xgenind@loc.fac), collapse="\n")

				# reformat allele data
					xtab=as.data.frame(xgenind@tab)

					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
						newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallels)
						
						if(msat==levels(xgenind@loc.fac)[1])
						{
							msatdfout=data.frame(newallels)
						} else
						{
							msatdfout=cbind(msatdfout, newallels)
						}
					}
					names(msatdfout)=levels(xgenind@loc.fac)

				# for each population
					popname=as.character(pop(xgenind))
					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
					popnameallels=paste(popname, allelsi, sep=" , ")
					
					outputlines=rbind(titlerow, locusnamesrow)
					for(popi in levels(pop(xgenind)))
					{
						outputlines=rbind(outputlines, "Pop")
						for(rowi in which(substr(popnameallels, 1,6)==popi))
						{
							outputlines=rbind(outputlines, popnameallels[rowi])
						}
					}
					
				return(outputlines)
					
			}
			print(" --> genind2genepop4diveRsity(genind-object, header, popstrata)")
# ... ... STRUCTURE
			genind2STRUCTURE <- function(geninddf,outputdir, datasetname, ORDERED=FALSE)
			{
				# adapt names and save data.names
					namesout=row.names(geninddf$tab)
					namesadapted=paste0(as.character(pop(setPop(geninddf,~Site))),"_",as.character(pop(setPop(geninddf,~AgeClass))))
					namesadapted=make.unique(namesadapted,sep="_")
					write.csv2(data.frame(Originalname=namesout, Nameadapted=namesadapted, Ageclass=as.character(pop(setPop(geninddf,~AgeClass))), HeightClass=as.character(pop(setPop(geninddf,~HeightClass)))),paste0(outputdir,datasetname,"_STRUCTURE_Names.csv"), row.names=FALSE, quote=FALSE)
					
				# convert allele counts to columns 
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(geninddf@loc.fac))
					{

						newallels=apply(geninddf$tab[,geninddf@loc.fac==msat],1,getallels)
						
						# split alleles to 2 cols
						# exchange "000" by "-9"
						newallels=gsub("000000", " -9 -9",newallels)
						
						if(msat==levels(geninddf@loc.fac)[1])
						{
							msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
						} else
						{
							msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
						}
					}
					names(msatdfout)=levels(geninddf@loc.fac)
					
					# add names
					msatdfout$Name=namesadapted
					
					# re-order
					msatdfoutsorted=msatdfout[,c(dim(msatdfout)[2],1:(dim(msatdfout)[2]-1))]
			
				# order individuals by populations
						if(ORDERED)
						{
							allpopi=as.character(pop(setPop(geninddf,~Site)))
							
							orderi=NULL
							for(popi in allpopsordered[which(allpopsordered%in%unique(allpopi))])
							{
								orderi=c(orderi,which(allpopi==popi))
							}
							msatdfoutsorted=msatdfoutsorted[orderi,]
						}
						
				# output
					write.table(msatdfoutsorted, paste0(outputdir,datasetname,"_STRUCTURE_Data.txt"), row.names=F, quote=F,col.names=F)

			}
			print(" --> genind2STRUCTURE(genind-object,outputdir,datasetname, ORDERED)")
			genind2STRUCTUREmitlocdat <- function(geninddf,outputdir, datasetname)
			{
				# adapt names and save data.names
					namesout=row.names(geninddf$tab)
					namesadapted=paste0(as.character(pop(setPop(geninddf,~Site))),"_",as.character(pop(setPop(geninddf,~AgeClass))))
					namesadapted=make.unique(namesadapted,sep="_")
					popnums=as.numeric(pop(setPop(geninddf,~Site)))
					write.csv2(data.frame(Originalname=namesout, Originalpop=as.character(pop(setPop(geninddf,~Site))), PopIDinSTRUCTUREdf=popnums, Nameadapted=namesadapted, Ageclass=as.character(pop(setPop(geninddf,~AgeClass))), HeightClass=as.character(pop(setPop(geninddf,~HeightClass)))),paste0(outputdir,datasetname,"_STRUCTURE_Names.csv"), row.names=FALSE, quote=FALSE)

				# convert allel counts to allel colums 
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(geninddf@loc.fac))
					{

						newallels=apply(geninddf$tab[,geninddf@loc.fac==msat],1,getallels)
						
						# split alleles to 2 cols
						# exchange "000" by "-9"
						newallels=gsub("000000", " -9 -9",newallels)
						
						if(msat==levels(geninddf@loc.fac)[1])
						{
							msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
						} else
						{
							msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
						}
					}
					names(msatdfout)=levels(geninddf@loc.fac)
					
					# add names
					msatdfout$Name=namesadapted
					msatdfout$Pop=popnums
					# re-order
					msatdfoutsorted=msatdfout[,c(which(names(msatdfout)=="Name"), which(names(msatdfout)=="Pop"),1:(dim(msatdfout)[2]-2))]
			
				# write output
					write.table(msatdfoutsorted, paste0(outputdir,datasetname,"_STRUCTURE_Data.txt"), row.names=F, quote=F,col.names=F)

			}
			print(" --> genind2STRUCTUREmitlocdat(genind-object,outputdir,datasetname)")
			genind2STRUCTURE2rows <- function(geninddf,outputdir, genindname)
			{
				# sort by population and age
					msattabsort=geninddf$tab[order(as.character(pop(setPop(geninddf,~Site/AgeClass)))),]
					popiclusiveageclass=as.character(pop(setPop(geninddf,~Site/AgeClass)))[order(as.character(pop(setPop(geninddf,~Site/AgeClass))))]
					popexclusiveageclass=as.character(pop(setPop(geninddf,~Site)))[order(as.character(pop(setPop(geninddf,~Site/AgeClass))))]
					
				# save strata
					write.csv2(data.frame(INDI=row.names(msattabsort), POP=popexclusiveageclass, POPAGE=popiclusiveageclass), paste0(outputdir,genindname,"_STRUCTURE_pops_ORIGINAL.csv"), row.names=FALSE, quote=FALSE)

				# convert allel counts to allel colums 
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(geninddf@loc.fac))
					{

						newallels=apply(msattabsort[,geninddf@loc.fac==msat],1,getallels)
						
						# split alleles to 2 cols
						# exchange "000" by "-9"
						newallels=gsub("000000", "-9 -9 ",newallels)
						
						if(msat==levels(geninddf@loc.fac)[1])
						{
							msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
						} else
						{
							msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
						}
					}
					namesout=NULL
					for(lvli in levels(geninddf@loc.fac))
					{
						namesout=c(namesout,paste0(lvli,c(",1",",2")))
					}
					names(msatdfout)=namesout
					
					# convert to 2-row data for each individual
						msatdfout$indinr=1:dim(msatdfout)[1]
						x1=msatdfout[,c(paste0(levels(geninddf@loc.fac),",1"),"indinr")]
						x2=msatdfout[,c(paste0(levels(geninddf@loc.fac),",2"),"indinr")]
						names(x2)=names(x1)
						msatoutmerged=rbind(x1,x2)
						msatoutmerged=msatoutmerged[order(msatoutmerged$indinr),]
			
				# add population
						iddf=rbind(data.frame(id=1:dim(msatdfout)[1], PopulationID=popexclusiveageclass),data.frame(id=1:dim(msatdfout)[1], PopulationID=popexclusiveageclass))
						msatoutmerged$PopulationID=iddf[order(iddf$id),]$PopulationID
						names(msatoutmerged)[which(names(msatoutmerged)=="indinr")]="IndividualID"
						msatoutmerged=msatoutmerged[,c((dim(msatoutmerged)[2]-1):(dim(msatoutmerged)[2]), 1:(dim(msatoutmerged)[2]-2))]
						msatoutmerged$IndividualID=paste0("Ind",msatoutmerged$IndividualID)
						names(msatoutmerged)=gsub(",1", "", names(msatoutmerged))
						
				# output
					write.table(msatoutmerged, paste0(outputdir,genindname,"_STRUCTURE2rows_individuals.txt"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
							
			}
			print(" --> genind2STRUCTURE2rows(genind-object,outputdir,genindname)")
# ... ... EEMS
			genind2EEMS <- function(geninddf,pathtofiles)
			{
				# adapt names and save data.names
					namesout=row.names(geninddf$tab)
					popsout=as.character(pop(setPop(geninddf,~Site)))
					
				#convert allele counts to allele columns 
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(geninddf@loc.fac))
					{
						newallels=apply(geninddf$tab[,geninddf@loc.fac==msat],1,getallels)
						
						# split alleles to 2 cols
						# exchange "000" by "-9"9
						newallels=gsub("000000", "-99-99",newallels)
						
						if(msat==levels(geninddf@loc.fac)[1])
						{
							msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
						} else
						{
							msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
						}
					}
					lvlout=NULL;for(i in 1:length(levels(geninddf@loc.fac)))lvlout=c(lvlout,rep(levels(geninddf@loc.fac)[i],2))
					names(msatdfout)=lvlout
					
					# add names
					dfout=cbind(Name=namesout, Population=popsout, msatdfout, geninddf$strata, other(geninddf))
					
					
				# check coordinates
						coord=as.data.frame(other(geninddf))[,c(2,1)]
							# ... switch if lat>90 
							poslat90=which(coord[,2]>90)
							cordbuffer=coord
							coord[poslat90,1]=cordbuffer[poslat90,2]
							coord[poslat90,2]=cordbuffer[poslat90,1]
							print(paste0("Position error, latitude switched with longitude: ", paste0(poslat90,collapse="-")))
							
				# define outer boundary as +1° of min/max coordinates
					londif=1
					latdif=1
					outer=data.frame(Lon=c(min(coord[,1])-londif,min(coord[,1])-londif,max(coord[,1])+londif,max(coord[,1])+londif,min(coord[,1])-londif), Lat=c(max(coord[,2])+latdif,min(coord[,2])-latdif,min(coord[,2])-latdif,max(coord[,2])+latdif,max(coord[,2])+latdif))
					
				# plot outer boundary and coordinates
					png(paste0(pathtofiles,"coord_outer.png"))
						plot(coord, pch=21, col="gray30", bg="gray70", xlim=c(min(coord[,1])-2*londif, max(coord[,1])+2*londif), ylim=c(min(coord[,2])-2*latdif, max(coord[,2])+2*latdif))
						polygon(outer,lty=2)
					dev.off()
						
				# output
					return(list(sites=msatdfout,coord, outer, dfout))

			}
			print(" --> genind2EEMS(genind-object,pathtofiles)")
# ... ... SPAGeDi
			genind2spagedi <- function(xgenind)
			{

				titlerow=paste(nInd(xgenind), 1, 2, nLoc(xgenind), 3, 2, sep="\t")
				locusnamesrow=paste(levels(xgenind@loc.fac), collapse="\n")

				# convert alleles 
					xtab=as.data.frame(xgenind@tab)
					
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
						newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallels)
						
						if(msat==levels(xgenind@loc.fac)[1])
						{
							msatdfout=data.frame(newallels)
						} else
						{
							msatdfout=cbind(msatdfout, newallels)
						}
					}
					names(msatdfout)=levels(xgenind@loc.fac)

				# for each population
					popname=as.character(pop(xgenind))
					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
					popnameallels=paste(popname, allelsi, sep=" , ")
					
					outputlines=rbind(titlerow, locusnamesrow)
					for(popi in levels(pop(xgenind)))
					{
						outputlines=rbind(outputlines, "Pop")
						for(rowi in which(substr(popnameallels, 1,6)==popi))
						{
							outputlines=rbind(outputlines, popnameallels[rowi])
						}
					}
					
					return(outputlines)
					
			}
			print(" --> genind2spagedi(genind-object)")
# ... ... TESS2
			genind2TESS2geno <- function(geninddf,outputdir, genindname, rmnocoords=FALSE)
			{
				# adapt names and save data.names
					namesoutindi=row.names(geninddf$tab)
					namesadapted=paste0(as.character(pop(setPop(geninddf,~Site))),"_",as.character(pop(setPop(geninddf,~AgeClass))))
					namesadapted=make.unique(namesadapted,sep="_")
					
				# convert allele counts to allele columns 
					getallels <- function(genindallels)
					{
						if(is.na(sum(genindallels))==TRUE)
						{
							return("000000")
						} else
						{
							if(sum(genindallels)==0)
							{
								return("000000")
							} else
							{
								allelselect=substr(names(genindallels),6,8)[which(genindallels!=0)]
								if(length(allelselect)==1)
								{
									return(paste0(allelselect,allelselect))
								} else
								{
									return(paste0(allelselect[1],allelselect[2]))
								}
							}
						}
					}
						
					# convert alleles to 6-characters, e.g. "000000"
						msatdfout=NULL
						for(msat in levels(geninddf@loc.fac))
						{

							newallels=apply(geninddf$tab[,geninddf@loc.fac==msat],1,getallels)
							
							# split alleles to 2 cols
							# exchange "000" by "-9"
							newallels=gsub("000000", "-9 -9 ",newallels)
							
							if(msat==levels(geninddf@loc.fac)[1])
							{
								msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
							} else
							{
								msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
							}
						}
						namesout=NULL
						for(lvli in levels(geninddf@loc.fac))
						{
							namesout=c(namesout,paste0(lvli,c(",1",",2")))
						}
						names(msatdfout)=namesout
						
				# add coordinates
					xydf=data.frame(X=gsub("NA", "00.000000", paste(as.vector(geninddf$other$latlong[,2]))),Y=gsub("NA", "00.000000", paste(as.vector(geninddf$other$latlong[,1]))))
					
				# re-order
					msatdfoutsorted=cbind(data.frame(Name=namesadapted), xydf, msatdfout)

				# output
					if(rmnocoords)
					{
						dfouti=msatdfoutsorted[order(pop(geninddf)),]
						write.table(dfouti[dfouti$X!="00.000000",], paste0(outputdir,genindname,".geno"), row.names=F, quote=F, col.names=TRUE)
						
						namout=data.frame(Originalname=namesoutindi, Nameadapted=namesadapted)[order(pop(geninddf)),]
						write.csv2(namout[dfouti$X!="00.000000",],paste0(outputdir,genindname,"_TESS2_Names.csv"), row.names=FALSE, quote=FALSE)

					} else {
						write.table(msatdfoutsorted[order(pop(geninddf)),], paste0(outputdir,genindname,".geno"), row.names=F, quote=F, col.names=TRUE)
						write.csv2(data.frame(Originalname=namesoutindi, Nameadapted=namesadapted)[order(pop(geninddf)),],paste0(outputdir,genindname,"_TESS2_Names.csv"), row.names=FALSE, quote=FALSE)
					}
			}
			print(" --> genind2TESS2geno(genind-object, outputdir, genindname, rmnocoords)")
# ... ... TESS3
			genind2TESS3geno <- function(genindf, genindname)
			{
				# access allele data
					coordsi=gsub("NA", "00.000000", paste(as.vector(unlist(apply(genindf$other$latlong[,2:1],1, function(x)paste(x,collapse=" "))))))[order(pop(genindf))]
					
					allelsi=gsub("NA", "9", as.vector(unlist(apply(as.data.frame(genindf$tab), 1, function(x)paste(x,collapse="")))))[order(pop(genindf))]
					
				# output
					write.table(coordsi, paste0(arbeitsdir, genindname,".coord"), row.names=FALSE, col.names=FALSE, quote=FALSE)
					write.table(allelsi, paste0(arbeitsdir, genindname,".geno"), row.names=FALSE, col.names=FALSE, quote=FALSE)
			}
			print(" --> genind2TESS3geno(genind-object,genindname)")
# ... ... BAPS
			genind2BAPS <- function(geninddf,outputdir,genindname)
			{
					
				# order by population and age
					msattabsort=geninddf$tab[order(as.character(pop(setPop(geninddf,~Site/AgeClass)))),]
					popiclusiveageclass=as.character(pop(setPop(geninddf,~Site/AgeClass)))[order(as.character(pop(setPop(geninddf,~Site/AgeClass))))]
					popexclusiveageclass=as.character(pop(setPop(geninddf,~Site)))[order(as.character(pop(setPop(geninddf,~Site/AgeClass))))]
					
					# reformat population identifier
						popindex=NULL
						for(popi in unique(popexclusiveageclass))
						{
							popindex=c(popindex, min(which(popexclusiveageclass==popi)))
						
						}
						write.table(data.frame(popindex=popindex), paste0(outputdir,genindname,"_BAPS_popindex.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
						
					# population name output
						write.table(data.frame(popnames=unique(popexclusiveageclass)), paste0(outputdir,genindname,"_BAPS_popnames.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
						
						# original population identifier output
							write.csv2(data.frame(INDI=row.names(msattabsort), POP=popexclusiveageclass, POPAGE=popiclusiveageclass), paste0(outputdir,genindname,"_BAPS_pops_ORIGINAL.txt"))

				# get coordinates and save them
					xydf=data.frame(X=gsub("NA", "0", paste(as.vector(geninddf$other$latlong[,2]))),Y=gsub("NA", "0", paste(as.vector(geninddf$other$latlong[,1]))))
					write.table(xydf, paste0(outputdir,genindname,"_BAPS_coords.txt"), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

				#convert allele counts to allele columns 
					# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(geninddf@loc.fac))
					{
						newallels=apply(msattabsort[,geninddf@loc.fac==msat],1,getallels)
						
						# split alleles to 2 cols
						# exchange "000" by "-9"
						newallels=gsub("000000", "-9 -9 ",newallels)
						
						if(msat==levels(geninddf@loc.fac)[1])
						{
							msatdfout=data.frame(cbind(substr(newallels,1,3),substr(newallels,4,6)))
						} else
						{
							msatdfout=cbind(msatdfout, cbind(substr(newallels,1,3),substr(newallels,4,6)))
						}
					}
					namesout=NULL
					for(lvli in levels(geninddf@loc.fac))
					{
						namesout=c(namesout,paste0(lvli,c(",1",",2")))
					}
					names(msatdfout)=namesout
					
					# convert to 2-row data
						msatdfout$indinr=1:dim(msatdfout)[1]
						x1=msatdfout[,c(paste0(levels(geninddf@loc.fac),",1"),"indinr")]
						x2=msatdfout[,c(paste0(levels(geninddf@loc.fac),",2"),"indinr")]
						names(x2)=names(x1)
						msatoutmerged=rbind(x1,x2)
						msatoutmerged=msatoutmerged[order(msatoutmerged$indinr),]
			
				# output
					write.table(msatoutmerged, paste0(outputdir,genindname,"_BAPS_individuals.txt"), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
							
			}
			print(" --> genind2BAPS(genind-object,outputdir,genindname)")
# ... ... CERVUS
			genind2genepop4CERVUS <- function(xgenind, header)
			{

				titlerow=header
				locusnamesrow=paste(levels(xgenind@loc.fac), collapse="\n")

				# adapt names and save data.names
					namesout=row.names(xgenind$tab)
					namesadapted=paste0(as.character(pop(setPop(xgenind,~Site))),"_",as.character(pop(setPop(xgenind,~AgeClass))))
					namesadapted=make.unique(namesadapted,sep="_")
					namesidsdf=data.frame(Originalname=namesout, Nameadapted=namesadapted, Ageclass=as.character(pop(setPop(xgenind,~AgeClass))), HeightClass=as.character(pop(setPop(xgenind,~HeightClass))))

				#convert allele counts to allele columns 
					xtab=as.data.frame(xgenind@tab)
					# convert alleles to 6-characters, e.g. "000000"
						msatdfout=NULL
						for(msat in levels(xgenind@loc.fac))
						{
								newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallels)
								
								if(msat==levels(xgenind@loc.fac)[1])
								{
										msatdfout=data.frame(newallels)
								} else
								{
										msatdfout=cbind(msatdfout, newallels)
								}
						}
						names(msatdfout)=levels(xgenind@loc.fac)

				# for each population
					popname=as.character(pop(xgenind))
					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
					popnameallels=paste(popname, allelsi, sep=" , ")
					
					outputlines=rbind(titlerow, locusnamesrow)
					for(popi in levels(pop(xgenind)))
					{
							outputlines=rbind(outputlines, "Pop")
							for(rowi in which(popname==popi))
							{
									outputlines=rbind(outputlines, popnameallels[rowi])
							}
					}
					
				return(list(outputlines,namesidsdf))

			}
			print(" --> genind2genepop4CERVUS(genind-object,header)")
			getallelsCervus <- function(genindallels)
			{
				if(is.na(sum(genindallels))==TRUE)
				{
					return("0,0")
				} else
				{
					if(sum(genindallels)==0)
					{
						return("0,0")
					} else
					{
						allelselect=substr(names(genindallels),6,8)[which(genindallels!=0)]
						if(length(allelselect)==1)
						{
							return(paste0(allelselect,",",allelselect))
						} else
						{
							return(paste0(allelselect[1],",",allelselect[2]))
						}
					}
				}
			}
			print(" --> getallelsCervus(genindallels)")
			genind2CERVUS_parentage <- function(xgenind)
			{
				locis=NULL;for(lvli in levels(xgenind@loc.fac))locis=c(locis, paste0(lvli,c("a","b")))
				locusnamesrow=paste(c("Individual ID", locis), collapse=",")
					
				# adapt names and save data.names
					namesout=row.names(xgenind$tab)
					namesadapted=paste0(as.character(pop(setPop(xgenind,~Site))),"_",as.character(pop(setPop(xgenind,~AgeClass))))
					namesadapted=make.unique(namesadapted,sep="_")
					namesidsdf=data.frame(Originalname=namesout, Nameadapted=namesadapted, Ageclass=as.character(pop(setPop(xgenind,~AgeClass))), HeightClass=as.character(pop(setPop(xgenind,~HeightClass))), Height=as.numeric(as.character(pop(setPop(xgenind,~Height)))), Site=as.character(pop(setPop(xgenind,~Site))))

				# genotype file
					xtab=as.data.frame(xgenind@tab)
						# nehme jeden Msat und mache 6stellige Angabe der beiden Allele wie z.B. "0" draus
							msatdfout=NULL
							for(msat in levels(xgenind@loc.fac))
							{
									# print(msat)
									
									newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallelsCervus)
									
									if(msat==levels(xgenind@loc.fac)[1])
									{#neu
											msatdfout=data.frame(newallels)
									} else
									{
											msatdfout=cbind(msatdfout, newallels)
									}
							}
							names(msatdfout)=levels(xgenind@loc.fac)

					# for each population
						allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=",")))
						popnameallels=paste(namesadapted, allelsi, sep=",")
						genotypefile=c(locusnamesrow,popnameallels)	
				
				# offspring file 
						headeroff=c("Offspring ID")
						offspringids=as.character(namesidsdf[namesidsdf$HeightClass=="Seed",]$Nameadapted)
						offspringfile=c(headeroff, offspringids)

						offspringids=as.character(namesidsdf[namesidsdf$HeightClass=="Seed" | namesidsdf$HeightClass=="Sapl",]$Nameadapted)
						offspringfileseedandsapl=c(headeroff, offspringids)
				# candidate parent file
					# .... one column for all offspring
						parentdfile=as.character(namesidsdf[namesidsdf$HeightClass!="Seed",]$Nameadapted)
						parentdfileheights=as.character(namesidsdf[namesidsdf$Height>=200 & !namesidsdf$Site%in%c("LLR","LLL","KO05","KO022"),]$Nameadapted)

				return(list(genotypefile,offspringfile,parentdfile,namesidsdf, parentdfileheights, offspringfileseedandsapl))
				
			}
			print(" --> genind2CERVUS_parentage(genind-object)")
# ... ... Geneland
			# ... internal function to converse allele-tab to merged allele sizes per locus
			getallelswithspace <- function(genindallels)
			{
				if(is.na(sum(genindallels))==TRUE)
				{
					return("000 000")
				} else
				{
					if(sum(genindallels)==0)
					{
						return("000 000")
					} else
					{
						allelselect=substr(names(genindallels),6,8)[which(genindallels!=0)]
						if(length(allelselect)==1)
						{
							return(paste0(allelselect," ",allelselect))
						} else
						{
							return(paste0(allelselect[1]," ",allelselect[2]))
						}
					}
				}
			}
			print(" --> getallelswithspace(genindallels)")
			genind2geneland <- function(xgenind)
			{
				xtab=as.data.frame(xgenind@tab)
					
				# convert alleles to 6-characters, e.g. "000000"
					msatdfout=NULL
					for(msat in levels(xgenind@loc.fac))
					{
							newallels=apply(xtab[,xgenind@loc.fac==msat],1,getallelswithspace)
							
							if(msat==levels(xgenind@loc.fac)[1])
							{
									msatdfout=data.frame(newallels)
							} else
							{
									msatdfout=cbind(msatdfout, newallels)
							}
					}
					names(msatdfout)=levels(xgenind@loc.fac)

					allelsi=as.vector(apply(msatdfout, 1, function(x)paste(x,collapse=" ")))
				
				return(allelsi)
							
			}
			print(" --> genind2geneland(genind-object)")

			


print(" --- loading finished! --- ")














