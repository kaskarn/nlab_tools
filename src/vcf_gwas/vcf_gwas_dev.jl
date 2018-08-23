using CodecZlib, GeneticVariation, DataFrames, CSV, GLM
#using VCFTools
#incmd = "--phenid famid --rhs age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+sitej --phen qrs --vcfpath /proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS/ARIC_AA_gwas_frz3/ARIC_AA_gwas_frz3.chr16.vcf.gz --phepath /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/phenotypes/ARIC_1kgp3/aric_ecg_tosugen_AA.txt --outfile /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/runs/test_julia/testout.txt"
#incmd = "--phenid analysis_id --rhs age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10 --phen ldl --vcfpath /proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS/ARIC_AA_gwas_frz3/ARIC_AA_gwas_frz3.chr16.vcf.gz --phepath /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/phenotypes/WHI_1kgp3/page_phen/pageii_whi_pheno_qt_qrs_ldl.txt --outfile /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/runs/test_julia/testout.txt"
incmd = join(ARGS, " ")
pinputs = Dict()
[get!(pinputs, split(i)[1], size(split(i),1) > 1 ? split(i)[2] : true) for i in split(incmd, "--")[2:end]]
  
[println("$(k): $(pinputs[k])") for k in keys(pinputs)];

#reader=VCF.Reader(GzipDecompressorStream(open(pinputs["vcfpath"])))
reader = VCFTools.VCF.Reader(openvcf(pinputs["vcfpath"]))
nsamp = nsamples(pinputs["vcfpath"])
g=NullableArray(zeros(nsamples)
@time for i in 1:10000
	copy_gt!(g, reader; model = :additive, impute = false, center = false, scale = false)
end
gw=NullableArray(zeros(nsamp,10000))
@time copy_gt!(gw, reader; model = :additive, impute = false, center = false, scale = false)

pheno = CSV.read(pinputs["phepath"], delim = '\t') #WAITING FOR UPDATE TO PARSER...

vcfids = DataFrame(id_rdy = reader.header.sampleID, ind = 1:length(reader.header.sampleID))
rename!(pheno, Symbol(pinputs["phenid"]) => :id_rdy)
vcf_phen = join(pheno, vcfids, on = :id_rdy, kind = :inner)
sort!(vcf_phen , cols = order(:ind))
vcfind = vcf_phen[:ind]

formula = Formula(Symbol(pinputs["phen"]), parse(join(["G + ", pinputs["rhs"]])))

outres = open(pinputs["outfile"], "w+")
write(outres, "VCF_ID\tCHROM\tPOS\tREF\tALT\tALT_AF\tALT_AC\tN_INFORMATIVE\tBETA\tSE\tPVALUE")

gind = trues(length(vcfind))
vcf_phen[:G] = rand(size(vcf_phen,1))
tmdf = fit(LinearModel, formula, vcf_phen).mf.df
mf = convert(Matrix{Float64}, Array(delete!(tmdf,1)))
y = vcf_phen[Symbol(pinputs["phen"])] 

function process_var!(vcfnow, y, mf, gind, reader, vcfind, outres)
	gnow = tryparse.(Float64, VCF.genotype(vcfnow, vcfind, "DS"))
	n, ac = 0, 0.0
	for (i, g) in enumerate(gnow)
		gind[i] = g.hasvalue
		if g.hasvalue
			n += 1
			ac += g.value
			mf[i,1] = g.value
		end
	end
	
	write(outres, '\n')
	(ac < 2.0 || ac > 2n-2) && return join(outres, [VCF.id(vcfnow)[1], VCF.chrom(vcfnow), VCF.pos(vcfnow), VCF.ref(vcfnow), VCF.alt(vcfnow)[1], ac/2/n, ac, n, "NA", "NA", "NA"], '\t')
	mnow = lm(mf[gind,:],y[gind])
	b = coef(mnow)[1]
	se = stderr(mnow)[1]
	p = GLM.ccdf(GLM.FDist(1,GLM.df_residual(mnow)),abs2(b/se))
	return join(outres, [VCF.id(vcfnow)[1], VCF.chrom(vcfnow), VCF.pos(vcfnow), VCF.ref(vcfnow), VCF.alt(vcfnow)[1], ac/2/n, ac, n, b, se, p], '\t')
end

println("GWAS Started on $(size(vcf_phen, 1)) Individuals\n\n")

i = 0
haskey(pinputs, "test") && (pinputs["test"] = parse(Int64, pinputs["test"]))
!haskey(pinputs, "test") && (get!(pinputs, "test", typemax(Int64)))
@time for vcfnow in Iterators.take(reader, pinputs["test"])
	i += 1
	process_var!(vcfnow, y, mf, gind, reader, vcfind, outres)
	i%10000 == 0 && println("$(i) variants processed")
end
