#####################################compare_genomics.script#####################################
01.准备数据
###建立工作目录
mkdir 01.prepare_data
cd 01.prepare_data
###软件安装###
###mamba install -c bioconda agat###
species=species
genome=species.fa
gff=species.gff.gz
#保留蛋白编码基因
agat_sp_filter_feature_by_attribute_value.pl --gff  $gff --attribute gene_biotype --value protein_coding -t '!' -o $species.protein_coding.gff3
#保留最长转录本
agat_sp_keep_longest_isoform.pl --gff $species.protein_coding.gff3 -o $species.longest_isoform.gff3
#提取cds序列
/home/liuxuan/miniconda3/envs/liuxuan-genome-funannotate/bin/gff3_file_to_proteins.pl  --gff3  $species.longest_isoform.gff3 --fasta $genome  --seqType CDS >$species.cds.fa
#提取pep序列
/home/liuxuan/miniconda3/envs/liuxuan-genome-funannotate/bin/gff3_file_to_proteins.pl  --gff3  $species.longest_isoform.gff3  --fasta  $genome --seqType prot >$species.pep.fa
#提取cds 转换成OrthoFinder 需要的序列格式
perl get_gene_longest_fa_for_OrthoFinder.pl -l 150 --gff $species.longest_isoform.gff3 --fa $species.cds.fa -p $species -o cds
#提取pep 转换成OrthoFinder 需要的序列格式
perl get_gene_longest_fa_for_OrthoFinder.pl -l 50 --gff $species.longest_isoform.gff3 --fa $species.pep.fa -p $species -o pep

02.基因家族聚类分析-同源基因鉴定
###建立工作目录
mkdir 02.gene_family
cd 02.gene_family
cp -r  ../01.prepare_data/pep ./
###软件安装###
#mamba install -c bioconda orthofinder
###运行orthofinder### >5h
nohup orthofinder  -f pep  -S diamond &

03.基因家族聚类分析-进化树的构建
###建立工作目录
mkdir 03.Phylo_Tree
cd 03.Phylo_Tree
###数据准备###
orthDir=../02.gene_family/pep/OrthoFinder/Results_*
cp   $orthDir/Orthogroups/Orthogroups.tsv ./
cp   $orthDir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ./
sed -i 's#\r##'  Orthogroups.tsv
###生成单拷贝基因家族列表文件###
###提取单拷贝基因家族成员###
awk '{if(NR==1){f=NF}; if( $0 !~ /,/ && NF==f ){print $0}}' Orthogroups.tsv > single_copy.txt
###对基因ID添加物种前缀,方便后续合并基因序列后用物种名字代表序列###
awk '{if(NR==1){n=split($0,A, "\t")} else { for(i=2; i<= n; i++){  printf A[i]"_"$i"\t" }; printf "\n" } }' single_copy.txt |sed 's/\s\+$//' > single_copy.txt.renamed
###多序列比对###
mkdir pep
for f in `ls ../01.prepare_data/pep/*fasta`;do
  filename=`basename $f`
  species=${filename%*.fasta}
  sed "s/>/>${species}_/" $f >pep/$filename
done

mkdir cds
for f in `ls ../01.prepare_data/cds/*fasta`;do
  filename=`basename $f`
  species=${filename%*.fasta}
  sed "s/>/>${species}_/" $f >cds/$filename
done
###合并所有物种的序列到一起###
cat ./pep/*fasta >all.pep.fasta
cat ./cds/*fasta >all.cds.fasta
###软件安装###
#mamba install -c bioconda  muscle=5.1
###ParaAT 多序列比对 5min###
perl ParaAT.pl -h single_copy.txt.renamed -a all.pep.fasta -n all.cds.fasta -o aln_out -p $threads -format fasta -m muscle -verbose
###CDS比对结果构建进化树###
###所有比对好的cds序列去掉基因名字部分，保留物种名部分，注意物种名字格式:*_*
ls aln_out/*.cds_aln.fasta |while read a;do sed  -r 's/([A-Za-z]+_[A-Za-z]+)_.*/\1/g'  $a > $a.renamed;done
###合并成supergene###
seqkit concat aln_out/*.cds_aln.fasta.renamed  >  single_copy.cds_msa.fasta
###剪切掉不保守的部分###
trimal -in single_copy.cds_msa.fasta -out single_copy.cds_msa_trimed.fasta -automated1
###运行raxml-ng### >24h
raxml-ng --all --msa single_copy.cds_msa_trimed.fasta --msa-format FASTA \
  --model GTR+G+I --prefix raxml.cds --threads 40 --seed $RANDOM --bs-trees 100 \
  --outgroup Oryza_sativa --tree pars{25},rand{25}

04.物种分歧时间估计
###建立工作目录###
mkdir 04.time_tree
cd 04.time_tree
###格式转换###
for f in ../03.Phylo_Tree/aln_out/*.cds_aln.fasta.renamed ;do 
 trimal  -in $f -out $f.phy  -phylip_paml -automated1
done
cat ../03.Phylo_Tree/aln_out/*.cds_aln.fasta.renamed.phy >input.phy
###获得基因数量###
ndata=`ls ../03.Phylo_Tree/aln_out/*.cds_aln.fasta.renamed|wc -l`
###复制进化树###
cp ../03.Phylo_Tree/raxml.cds.raxml.support input.tre ./
###进化树去掉不必要枝长和boot值###
cat species.tre |sed -r 's/:[0-9\.]+//g' |sed -r 's/\)[0-9\.]+/)/g'>species_formated.tre
###手动编写进化树
#添加化石时间点：timetree http://timetree.org/
#设置进化树的root : https://itol.embl.de/
#添加化石校准点时间信息（格式是时间范围’>0.23<0.26’或者时间点‘@0.245’)，单位时百万年前100Ma；
#再在首行添加两个数字（物种数量和树的数量），空格隔开，可得到input.tre文件。
###生成mcmctree配置文件：mcmctree.ctl###
echo "
          seed = -1        *设置随机数作为seed，-1代表使用系统当前时间作为随机数
       seqfile = input.phy *输入多序列比对文件
      treefile = input.tre *带校准点（化石时间）的有根树文件
       outfile = out.txt   *输出文件
      mcmcfile = mcmc.txt  *输出的mcmc信息文件，可用Tracer软件查看
         
         ndata = $ndata   * 输入的多序列比对的数据区域的数量；多个数据phy格式合并
       seqtype = 0   * 设置多序列比对数据类型；0：核酸数据；1：密码子比对数据；2：氨基酸数据；
       usedata = 3   * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
                     * 是否利用多序列比对数据；
                     * 0: no data,不使用，不会进行likelihood估算，会快速得到mcmc树，但分歧时间不可用; 
                     * 1:seq like，使用多序列比对数据进行likelihood估算，正常进行mcmc; usedata=1时model无法选择；
                     * 2:use in.BV, 进行正常的approximation likelihood分析，不读取多序列比对数据，直接读取当前目录的in.BV文件，in.BV是由usedata = 3时生成的out.BV重命名得来；此外，由于程序BUG，当设置usedata = 2时，一定要在改行参数后加 *，否则程序报错 Error: file name empty..；
                     * 3：out.BV,程序利用多序列比对数据调用baseml/codeml命令对数据进行分析，生成out.BV文件。由于mcmctree调用baseml/codeml进行估算的参数设置可能不太好（特别时对蛋白序列进行估算时），推荐自己修改软件自动生成的baseml/codeml配置文件，然后再手动运行baseml/codeml命令，再整合其结果文件为out.BV文件。
 
         clock = 2       * 设置分子钟算法，1: global clock，表示所有分支进化速率一致; 2: independent rates，各分支的进化速率独立且进化速率的对数log(r)符合正态分布; 3，correlated rates方法，和方法2类似，但是log(r)的方差和时间t相关。
                         *       TipDate = 1 100  *当外部节点由取样时间时使用该参数进行设置，同时该参数也设置了时间单位。具体数据示例请见examples/TipData文件夹。
        RootAge = '<1.73'  * constraint on root age, used if no fossil for root.设置root节点的分歧时间，一般设置一个最大值。

         model = 4      * models for DNA:
                        *  0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85；*设置碱基替换模型；当设置usedata = 1时，model不能使用超过4的模型，所以usedata = 1时用model = 4；usedata不等于1时，用model = 7，即GTR模型；
                        * models for codons:
                        *  0:one 恒定速率模型, 1:b 中性模型, 2:2 or more dN/dS ratios for branches 选择模型。
                        * models for AAs or codon-translated AAs:
                        *  0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        *  6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
         alpha = 0.5    * alpha for gamma rates at sites；*核酸序列中不同位点，其进化速率不一致，其变异速率服从GAMMA分布。一般设置GAMMA分布的alpha值为0.5。若该参数值设置为0，则表示所有位点的进化速率一致。此外，当userdata = 2时，alpha、ncatG、alpha_gamma、kappa_gamma这4个GAMMA参数无效。因为userdata = 2时，不会利用到多序列比对的数据。
         ncatG = 5      * No. categories in discrete gamma；设置离散型GAMMA分布的categories值。

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)? 
                      * 设置是否移除不明确字符（N、？、W、R和Y等）或含以后gap的列后再进行数据分析： 0，不需要，但在序列两两比较的时候，还是会去除后进行比较；  

       BDparas = 1 1 0.1   * birth, death, sampling；*设置出生率、死亡率和取样比例。若输入有根树文件中的时间单位发生改变，则需要相应修改出生率和死亡率的值。例如，时间单位由100Myr变换为1Myr，则要设置成'.01 .01 0.1'。
   kappa_gamma = 6 2       * gamma prior for kappa；设置kappa（转换/颠换比率）的GAMMA分布参数。
   alpha_gamma = 1 1       * gamma prior for alpha；设置GAMMA形状参数alpha的GAMMA分布参数。只对usedata=1时起作用，其他2，3不起作用

   rgene_gamma = 2 20 1   * alpba bu a prior, gamma prior for rate for genes；物种替换速率可查文献；设置序列中所有位点平均[碱基/密码子/氨基酸]替换率的Dirichlet-GAMMA分布参数：alpha=2、beta=20、初始平均替换率为每100million年（取决于输入有根树文件中的时间单位）1个替换。若时间单位由100Myr变换为1Myr，则要设置成'2 2000 1'。总体上的平均进化速率为：2/20, 0.1 个替换 / 每100Myr，即每个位点每年的替换数为 1e-9。
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)；设置所有位点进化速率取对数后方差（sigma的平方）的Dirichlet-GAMMA分布参数：alpha=1、beta=10、初始方差值为1。当clock参数值为1时，表示使用全局的进化速率，各分枝的进化速率没有差异，即方差为0，该参数无效；当clock参数值为2时，若修改了时间单位，该参数值不需要改变；当clock参数值为3时，若修改了时间单位，该参数值需要改变。

      finetune = 1: .1 .1 .1 .1 .01 .1  * times, rates, mixing, paras, RateParas；冒号前的值设置是否自动进行finetune，一般设置成1，然程序自动进行优化分析；冒号后面设置各个参数的步进值：times, musigma2, rates, mixing, paras, FossilErr。由于有了自动设置，该参数不像以前版本那么重要了，可能以后会取消该参数。

         print = 1      *设置打印mcmc的取样信息：0，不打印mcmc结果；1，打印除了分支进化速率的其它信息（各内部节点分歧时间、平均进化速率、sigma2值）；2，打印所有信息。 
        burnin = 4000   *将前4000次迭代burnin后，再进行取样（即打印出该次迭代估算的结果信息，各内部节点分歧时间、平均进化速率、sigma2值和各分支进化速率等）。
      sampfreq = 100      *每100次迭代则取样一次
       nsample = 200000  *当取样次数达到该次数时，则取样结束，同时结束程序。

*** Note: Make your window wider (100 columns) when running this program.
" >mcmctree1.ctl
###运行mcmctree###
mcmctree mcmctree1.ctl >run1.o
mv out.BV in.BV
###近似似然法计算###
###生成mcmctree配置文件：mcmctree2.ctl###
echo "
          seed = -1        *设置随机数作为seed，-1代表使用系统当前时间作为随机数
       seqfile = input.phy *输入多序列比对文件
      treefile = input.tre *带校准点（化石时间）的有根树文件
       outfile = out.txt   *输出文件
      mcmcfile = mcmc.txt  *输出的mcmc信息文件，可用Tracer软件查看
         
         ndata = 1   * 输入的多序列比对的数据区域的数量；多个数据phy格式合并
       seqtype = 0   * 设置多序列比对数据类型；0：核酸数据；1：密码子比对数据；2：氨基酸数据；
       usedata = 2   * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
                     * 是否利用多序列比对数据；
                     * 0: no data,不使用，不会进行likelihood估算，会快速得到mcmc树，但分歧时间不可用; 
                     * 1:seq like，使用多序列比对数据进行likelihood估算，正常进行mcmc; usedata=1时model无法选择；
                     * 2:use in.BV, 进行正常的approximation likelihood分析，不读取多序列比对数据，直接读取当前目录的in.BV文件，in.BV是由usedata = 3时生成的out.BV重命名得来；此外，由于程序BUG，当设置usedata = 2时，一定要在改行参数后加 *，否则程序报错 Error: file name empty..；
                     * 3：out.BV,程序利用多序列比对数据调用baseml/codeml命令对数据进行分析，生成out.BV文件。由于mcmctree调用baseml/codeml进行估算的参数设置可能不太好（特别时对蛋白序列进行估算时），推荐自己修改软件自动生成的baseml/codeml配置文件，然后再手动运行baseml/codeml命令，再整合其结果文件为out.BV文件。
 
         clock = 2       * 设置分子钟算法，1: global clock，表示所有分支进化速率一致; 2: independent rates，各分支的进化速率独立且进化速率的对数log(r)符合正态分布; 3，correlated rates方法，和方法2类似，但是log(r)的方差和时间t相关。
                         *       TipDate = 1 100  *当外部节点由取样时间时使用该参数进行设置，同时该参数也设置了时间单位。具体数据示例请见examples/TipData文件夹。
        RootAge = '<2'  * constraint on root age, used if no fossil for root.设置root节点的分歧时间，一般设置一个最大值。

         model = 4      * models for DNA:
                        *  0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85；*设置碱基替换模型；当设置usedata = 1时，model不能使用超过4的模型，所以usedata = 1时用model = 4；usedata不等于1时，用model = 7，即GTR模型；
                        * models for codons:
                        *  0:one 恒定速率模型, 1:b 中性模型, 2:2 or more dN/dS ratios for branches 选择模型。
                        * models for AAs or codon-translated AAs:
                        *  0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        *  6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
         alpha = 0.5    * alpha for gamma rates at sites；*核酸序列中不同位点，其进化速率不一致，其变异速率服从GAMMA分布。一般设置GAMMA分布的alpha值为0.5。若该参数值设置为0，则表示所有位点的进化速率一致。此外，当userdata = 2时，alpha、ncatG、alpha_gamma、kappa_gamma这4个GAMMA参数无效。因为userdata = 2时，不会利用到多序列比对的数据。
         ncatG = 5      * No. categories in discrete gamma；设置离散型GAMMA分布的categories值。

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)? 
                      * 设置是否移除不明确字符（N、？、W、R和Y等）或含以后gap的列后再进行数据分析： 0，不需要，但在序列两两比较的时候，还是会去除后进行比较；  

       BDparas = 1 1 0.1   * birth, death, sampling；*设置出生率、死亡率和取样比例。若输入有根树文件中的时间单位发生改变，则需要相应修改出生率和死亡率的值。例如，时间单位由100Myr变换为1Myr，则要设置成'.01 .01 0.1'。
   kappa_gamma = 6 2       * gamma prior for kappa；设置kappa（转换/颠换比率）的GAMMA分布参数。
   alpha_gamma = 1 1       * gamma prior for alpha；设置GAMMA形状参数alpha的GAMMA分布参数。只对usedata=1时起作用，其他2，3不起作用

   rgene_gamma = 2 20 1   * au bu a prior, gamma prior for rate for genes；物种替换速率可查文献；设置序列中所有位点平均[碱基/密码子/氨基酸]替换率的Dirichlet-GAMMA分布参数：alpha=2、beta=20、初始平均替换率为每100million年（取决于输入有根树文件中的时间单位）1个替换。若时间单位由100Myr变换为1Myr，则要设置成'2 2000 1'。总体上的平均进化速率为：1 个替换 / 每100Myr，即每个位点每年的替换数为 1e-8。
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)；设置所有位点进化速率取对数后方差（sigma的平方）的Dirichlet-GAMMA分布参数：alpha=1、beta=10、初始方差值为1。当clock参数值为1时，表示使用全局的进化速率，各分枝的进化速率没有差异，即方差为0，该参数无效；当clock参数值为2时，若修改了时间单位，该参数值不需要改变；当clock参数值为3时，若修改了时间单位，该参数值需要改变。

      finetune = 1: .1 .1 .1 .1 .01 .1  * times, rates, mixing, paras, RateParas；冒号前的值设置是否自动进行finetune，一般设置成1，然程序自动进行优化分析；冒号后面设置各个参数的步进值：times, musigma2, rates, mixing, paras, FossilErr。由于有了自动设置，该参数不像以前版本那么重要了，可能以后会取消该参数。

         print = 1      *设置打印mcmc的取样信息：0，不打印mcmc结果；1，打印除了分支进化速率的其它信息（各内部节点分歧时间、平均进化速率、sigma2值）；2，打印所有信息。 
        burnin = 4000   *将前4000次迭代burnin后，再进行取样（即打印出该次迭代估算的结果信息，各内部节点分歧时间、平均进化速率、sigma2值和各分支进化速率等）。
      sampfreq = 100      *每100次迭代则取样一次
       nsample = 200000  *当取样次数达到该次数时，则取样结束，同时结束程序。

*** Note: Make your window wider (100 columns) when running this program.
" >mcmctree2.ctl
###运行mcmctree###
mcmctree mcmctree2.ctl

05.基因家族扩张和收缩
###建立工作目录###
mkdir 05.gene_family_cafe5
cd 05.gene_family_cafe5
###复制基因家族在不同物种中的数量表格###
cp  ../02.gene_family/pep/OrthoFinder/Results_*/Orthogroups/Orthogroups.GeneCount.tsv ./
###去除windows字符###
sed -i 's#\r##'  Orthogroups.GeneCount.tsv #去掉windows换行符
###修改第一列列名###
sed 's/[a-zA-Z0-9]\+$//' Orthogroups.GeneCount.tsv | awk '{print $1"\t"$0}'  |sed 's/Orthogroup/Desc/'  > Orthogroups.GeneCount.cafe.tsv
###选择树：二叉的（binary），有根的（rooted），超度量(时间树，ultrametric)的newick格式树
cp  ../04.time_tree_BV/FigTree.tre ./
###进化树转换成NWK格式，并把中括号及里面的内容删除###
cat FigTree.tre  | sed  's/\[[^]]\+\]//g'| awk -F "=" '/UTREE/{print $2} '  > tmp.tree.nwk
###乘以100，转换时间单位百万年###
sed  -e 's/:/\n:/g' -e 's/\([),]\)/\n\1/g'  tmp.tree.nwk  |awk '{if($1~/:$/){printf ":"100*$2} else {printf $0}}' |sed 's/\s\+//g' > cafe_formated.nwk
###先手动安装cafe5### 因为conda安装，安装不了其他脚本
#git clone https://github.com/hahnlab/CAFE5.git
#cd CAFE5
#./configure
#make
#过滤掉在不同物种中数量变化太大的基因
python3 /home/liuxuan/software/cafe/CAFE5/tutorial/clade_and_size_filter.py -i  Orthogroups.GeneCount.cafe.tsv -o  Orthogroups.GeneCount.cafe.filtered.tsv -s
###再用conda安装cafe5### 
#mamba seach cafe
#mamba install -c bioconda cafe
###运行 cafe5 ### 5min
cafe5  --infile  Orthogroups.GeneCount.cafe.filtered.tsv --tree cafe_formated.nwk --output_prefix  cafe_result  --cores 20     -p --pvalue 0.05 &>base.log
cafe5  --infile  Orthogroups.GeneCount.cafe.filtered.tsv --tree cafe_formated.nwk --output_prefix  cafe_result2  --cores 20    -k 2 -p --pvalue 0.05 &>k2.log
cafe5  --infile  Orthogroups.GeneCount.cafe.filtered.tsv --tree cafe_formated.nwk --output_prefix  cafe_result3  --cores 20    -k 3 -p --pvalue 0.05 &>k3.log
cafe5  --infile  Orthogroups.GeneCount.cafe.filtered.tsv --tree cafe_formated.nwk --output_prefix  cafe_result4  --cores 20    -k 4 -p --pvalue 0.05 &>k4.log
cafe5  --infile  Orthogroups.GeneCount.cafe.filtered.tsv --tree cafe_formated.nwk --output_prefix  cafe_result5  --cores 20    -k 5 -p --pvalue 0.05 &>k5.log
##选择最优结果
###k值的结果比较###
###查看k2，k3，k4, k5，等不同的结果文件Gamma_results.txt文件中的第一行信息，Model Gamma Final Likelihood (-lnL)值，挑选最大的为最优结果。
find ./ -name Gamma_results.txt |xargs grep "Model Gamma Final Likelihood"
###提取其中一个树，并删除枝上的基因家族信息,删除*号，删除节点家族基因数量###
awk 'NR==3 {print $NF}' cafe_result2/*_asr.tre |sed  -r  's/_[0-9]+//g ; s/\*//g '  > cafe_result2/cafe5_ready.nwk
###基因家族总数###
GFnum=`grep  "TREE " cafe_result2/Gamma_asr.tre|wc -l`
export QT_QPA_PLATFORM="offscreen"
python3 cafe5_plot.py --nwk cafe_result2/cafe5_ready.nwk --clade_results  cafe_result2/Gamma_clade_results.txt --gfn $GFnum -d cafe_result2
###提取扩张和收缩基因###
#提取Gamma_change.tab第n列为物种Talaromyces_barcinensis的扩张/收缩的orthogroupsID
cat cafe_result2/Gamma_change.tab |cut -f1,n|awk '{if($2>0){print $0}}' >Talaromyces_barcinensis.expanded 
cat cafe_result2/Gamma_change.tab |cut -f1,n|awk '{if($2<0){print $0}}' >Talaromyces_barcinensis.contracted  
#提取显著扩张或收缩的orthogroupsID
cat cafe_result2/Gamma_family_results.txt |grep "y"|cut -f1 >Talaromyces_barcinensis.p0.05.significant 
#提取显著扩张/收缩的物种的orthogroupsID
grep -f Talaromyces_barcinensis.p0.05.significant Talaromyces_barcinensis.expanded |cut -f1>Talaromyces_barcinensis.expanded.significant 
grep -f Talaromyces_barcinensis.p0.05.significant Talaromyces_barcinensis.contracted |cut -f1 >Talaromyces_barcinensis.contracted.significant 
#提取显著扩张/收缩的基因列表
grep -f Talaromyces_barcinensis.expanded.significant ../02.gene_family/pep/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt |sed "s/ /\n/g"|grep "^gene" |sort  |uniq >Talaromyces_barcinensis.expanded.significant.genes 
grep -f Talaromyces_barcinensis.contracted.significant ../02.gene_family/pep/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt |sed "s/ /\n/g"|grep "^gene" |sort |uniq >Talaromyces_barcinensis.contracted.significant.genes 
#提取显著扩张/收缩的基因序列
seqkit grep -f Talaromyces_barcinensis.expanded.significant.genes ../01.prepare_data/pep/Talaromyces_barcinensis.fasta >Talaromyces_barcinensis.expanded.significant.pep.fa 
seqkit grep -f Talaromyces_barcinensis.contracted.significant.genes ../01.prepare_data/pep/Talaromyces_barcinensis.fasta >Talaromyces_barcinensis.contracted.significant.pep.fa 

其他代码
###元素循环脚本###
#碳循环#
perl CCycdb.PL -situation assembly-based -wd ./ -m diamond -f fasta -s prot -tpm 1 -norm 0 -thread 32 -od ./Talaromyces_Ccycle-gene
#氮循环#
perl NCycProfiler.PL NCycProfiler.pl -d ./ -m diamond -f fasta -s prot -si sample-size -o Talaromyces_Ncycle-gene
#硫循环#
perl SCycDB_FunctionProfiler.PL NCycProfiler.pl -d ./ -m diamond -f fasta -s prot -si sample-size -o Talaromyces_Scycle-gene
#甲烷循环#
perl MCycDB_FunctionProfiler.PL NCycProfiler.pl -d ./ -m diamond -f fasta -s prot -si sample-size -o Talaromyces_Mcycle-gene

###CAZYmes注释###
nohup run_dbcan CAZyme_annotation --input_raw_data Talaromyces_barcinensis.fasta --output_dir out_auto2 --db_dir db --mode protein --e_value_threshold 1e-102 --coverage_threshold_dbcan 0.35 --e_value_threshold_dbcan 1e-15 --coverage_threshold_dbsub 0.35 --e_value_threshold_dbsub 1e-15 &



