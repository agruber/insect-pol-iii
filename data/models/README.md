wget -O RF00003.cm https://rfam.org/family/RF00003/cm
wget -O RF00004.cm https://rfam.org/family/RF00004/cm
wget -O RF00007.cm https://rfam.org/family/RF00007/cm
wget -O RF00009.cm https://rfam.org/family/RF00009/cm
wget -O RF00012.cm https://rfam.org/family/RF00012/cm
wget -O RF00015.cm https://rfam.org/family/RF00015/cm
wget -O RF00020.cm https://rfam.org/family/RF00020/cm
wget -O RF00026.cm https://rfam.org/family/RF00026/cm
wget -O RF00030.cm https://rfam.org/family/RF00030/cm
wget -O RF00066.cm https://rfam.org/family/RF00066/cm
wget -O RF00096.cm https://rfam.org/family/RF00096/cm
wget -O RF00548.cm https://rfam.org/family/RF00548/cm
wget -O RF00618.cm https://rfam.org/family/RF00618/cm
wget -O RF00619.cm https://rfam.org/family/RF00619/cm
wget -O RF01052.cm https://rfam.org/family/RF01052/cm

wget -O S1.tar "https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1005383.s010"
tar -xf S1.tar eukaryota.stk
mv eukaryota.stk tRNAsec.stk
# edit name in tRNAsec.stk to tRNAsec

cmbuild -F tRNAsec.cm tRNAsec.stk
cmcalibrate --cpu 8 tRNAsec.cm

cat RF0*.cm tRNAsec.cm > snRNAtRNA.cm
cmpress snRNAtRNA.cm
