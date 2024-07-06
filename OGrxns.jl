using DelimitedFiles
using StaticArrays
include("TXTL_Parameters.jl")

function ogrxnsetup(S)

enz_data = readdlm("enzyme_name_conc_nM_kcat_1_over_s_rxn_no.dat")
E = enz_data[2,:]./1e6
E_idx = findall(E.<2.4e-5) #find all enzymes that were not reported
E[E_idx] .= 5e-5 #Set all unknown enzymes to 50nM
 E .= 1

kcat = enz_data[3,:].*3600;
kcat = readdlm("config/KineticDict/kcat_1over_sec.dat").*3600;

Km = readdlm("config/KineticDict/Km_mM_matrix.dat");

TXTL_parameters = zeros(4,1)
TXTL_parameters[1] = .0675; #RNAP_concentration  #uM # 60-75nM (ACS SynBio Garamella 2016)
TXTL_parameters[2] = 20;  #max_transcription_rate  # >5 NT/s (ACS SynBio Garamella 2016)
TXTL_parameters[3] = 2.15; #RIBOSOME_concentration #uM #0.0016mM with 72% MaxActive (Underwood, Swartz, Puglisi 2005 Biotech Bioeng) & <0.0023mM (ACS SynBio Garamella 2016)
TXTL_parameters[4] = 1.5; #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)

TXTL = TXTLDictionary(TXTL_parameters)

VL = TXTL[:VL]
length_factor_translation = TXTL[:length_factor_translation]
L_tau_factor = TXTL[:L_tau_factor]
KL = TXTL[:KL]
VX = TXTL[:VX]
length_factor_transcription = TXTL[:length_factor_transcription]
X_tau_factor = TXTL[:X_tau_factor]
KX = TXTL[:KX]

    #Calculate TXTL rates
    TX = compute_transcription_rate(TXTL)
    TL = compute_translation_rate(S[152],TXTL)
    mRNA_degradation_rate = TXTL[:kdX]
    km_tx = TXTL[:km_tx]
    km_tx = 0.00001


ogrxn = zeros(eltype(S),200)

     #Redefine maltodextrin consumption rate and TXTL rates
  
     ogrxn[1] = kcat[1]*E[1] * S[95]/(Km[95,1]+S[95])
# ogrxn[1] = kcat[1]*E[1] * S[95]/(Km[95,1]+S[95]) * S[116]/(Km[116,1]+S[116])
ogrxn[2] = kcat[2]*E[2]
ogrxn[3] = kcat[3]*E[3]
ogrxn[4] = kcat[4]*E[4]
ogrxn[5] = kcat[5]*E[5] * S[96]/(Km[96,5]+S[96])
ogrxn[6] = kcat[6]*E[6]
ogrxn[7] = kcat[7]*E[7] * S[31]/(Km[31,7]+S[31]) * S[61]/(Km[61,7]+S[61])
ogrxn[8] = kcat[8]*E[8] * S[57]/(Km[57,8]+S[57])
ogrxn[9] = kcat[9]*E[9] * S[48]/(Km[48,9]+S[48])
ogrxn[10] = kcat[10]*E[10] * S[31]/(Km[31,10]+S[31]) * S[48]/(Km[48,10]+S[48])
ogrxn[11] = kcat[11]*E[11] * S[50]/(Km[50,11]+S[50])
ogrxn[12] = kcat[12]*E[12] * S[50]/(Km[50,12]+S[50])
ogrxn[13] = kcat[13]*E[13] * S[56]/(Km[56,13]+S[56])
ogrxn[14] = kcat[14]*E[14]
ogrxn[15] = kcat[15]*E[15] * S[56]/(Km[56,15]+S[56])
ogrxn[16] = kcat[16]*E[16] * S[107]/(Km[107,16]+S[107])
ogrxn[17] = kcat[17]*E[17] * S[67]/(Km[67,17]+S[67]) * S[106]/(Km[106,17]+S[106])
ogrxn[18] = kcat[18]*E[18] * S[56]/(Km[56,18]+S[56]) * S[104]/(Km[104,18]+S[104])
ogrxn[19] = kcat[19]*E[19] * S[105]/(Km[105,19]+S[105])
ogrxn[20] = kcat[20]*E[20] * S[18]/(Km[18,20]+S[18])
ogrxn[21] = kcat[21]*E[21] * S[6]/(Km[6,21]+S[6]) * S[31]/(Km[31,21]+S[31])
ogrxn[22] = kcat[22]*E[22] * S[6]/(Km[6,22]+S[6])
ogrxn[23] = kcat[23]*E[23]
ogrxn[24] = kcat[24]*E[24]
ogrxn[25] = kcat[25]*E[25] * S[113]/(Km[113,25]+S[113])
ogrxn[26] = kcat[26]*E[26] * S[18]/(Km[18,26]+S[18]) * S[113]/(Km[113,26]+S[113])
ogrxn[27] = kcat[27]*E[27] * S[31]/(Km[31,27]+S[31]) * S[110]/(Km[110,27]+S[110])
ogrxn[28] = kcat[28]*E[28] * S[113]/(Km[113,28]+S[113])
ogrxn[29] = kcat[29]*E[29] * S[104]/(Km[104,29]+S[104]) * S[122]/(Km[122,29]+S[122])
ogrxn[30] = kcat[30]*E[30] * S[31]/(Km[31,30]+S[31]) * S[122]/(Km[122,30]+S[122])
ogrxn[31] = kcat[31]*E[31] * S[57]/(Km[57,31]+S[57]) * S[106]/(Km[106,31]+S[106])
ogrxn[32] = kcat[32]*E[32] * S[107]/(Km[107,32]+S[107])
ogrxn[33] = kcat[33]*E[33]
ogrxn[34] = kcat[34]*E[34] * S[11]/(Km[11,34]+S[11]) * S[106]/(Km[106,34]+S[106])
ogrxn[35] = kcat[35]*E[35] * S[126]/(Km[126,35]+S[126])
ogrxn[36] = kcat[36]*E[36]
ogrxn[37] = kcat[37]*E[37] * S[125]/(Km[125,37]+S[125])
ogrxn[38] = kcat[38]*E[38] * S[126]/(Km[126,38]+S[126])
ogrxn[39] = kcat[39]*E[39] * S[56]/(Km[56,39]+S[56]) * S[127]/(Km[127,39]+S[127])
ogrxn[40] = kcat[40]*E[40] * S[46]/(Km[46,40]+S[46]) * S[48]/(Km[48,40]+S[48])
ogrxn[41] = kcat[41]*E[41] * S[125]/(Km[125,41]+S[125])
ogrxn[42] = kcat[42]*E[42] * S[56]/(Km[56,42]+S[56]) * S[127]/(Km[127,42]+S[127])
ogrxn[43] = kcat[43]*E[43] * S[46]/(Km[46,43]+S[46])
ogrxn[44] = kcat[44]*E[44] * S[48]/(Km[48,44]+S[48]) * S[56]/(Km[56,44]+S[56])
ogrxn[45] = kcat[45]*E[45] * S[11]/(Km[11,45]+S[11])
ogrxn[46] = kcat[46]*E[46]
ogrxn[47] = kcat[47]*E[47] * S[16]/(Km[16,47]+S[16]) * S[110]/(Km[110,47]+S[110])
ogrxn[48] = kcat[48]*E[48] * S[36]/(Km[36,48]+S[36])
ogrxn[49] = kcat[49]*E[49] * S[81]/(Km[81,49]+S[81])
ogrxn[50] = kcat[50]*E[50] * S[81]/(Km[81,50]+S[81]) * S[106]/(Km[106,50]+S[106])
ogrxn[51] = kcat[51]*E[51] * S[21]/(Km[21,51]+S[21]) * S[107]/(Km[107,51]+S[107])
ogrxn[52] = kcat[52]*E[52] * S[21]/(Km[21,52]+S[21]) * S[104]/(Km[104,52]+S[104])
ogrxn[53] = kcat[53]*E[53] * S[18]/(Km[18,53]+S[18])
ogrxn[54] = kcat[54]*E[54] * S[131]/(Km[131,54]+S[131])
ogrxn[55] = kcat[55]*E[55] * S[54]/(Km[54,55]+S[54])
ogrxn[56] = kcat[56]*E[56] * S[54]/(Km[54,56]+S[54])
ogrxn[57] = kcat[57]*E[57] * S[91]/(Km[91,57]+S[91])
ogrxn[58] = kcat[58]*E[58] * S[91]/(Km[91,58]+S[91]) * S[104]/(Km[104,58]+S[104])
ogrxn[59] = kcat[59]*E[59] * S[105]/(Km[105,59]+S[105]) * S[110]/(Km[110,59]+S[110])
ogrxn[60] = kcat[60]*E[60]
ogrxn[61] = kcat[61]*E[61]
ogrxn[62] = kcat[62]*E[62]
ogrxn[63] = kcat[63]*E[63] * S[18]/(Km[18,63]+S[18])
ogrxn[64] = kcat[64]*E[64]
ogrxn[65] = kcat[65]*E[65] * S[105]/(Km[105,65]+S[105])
ogrxn[66] = kcat[66]*E[66] * S[104]/(Km[104,66]+S[104]) * S[107]/(Km[107,66]+S[107])
ogrxn[67] = kcat[67]*E[67] * S[105]/(Km[105,67]+S[105]) * S[106]/(Km[106,67]+S[106])
ogrxn[68] = kcat[68]*E[68] * S[105]/(Km[105,68]+S[105])
ogrxn[69] = kcat[69]*E[69] * S[105]/(Km[105,69]+S[105])
ogrxn[70] = kcat[70]*E[70]
ogrxn[71] = kcat[71]*E[71] * S[81]/(Km[81,71]+S[81])
ogrxn[72] = kcat[72]*E[72] * S[16]/(Km[16,72]+S[16])
ogrxn[73] = kcat[73]*E[73] * S[91]/(Km[91,73]+S[91]) * S[104]/(Km[104,73]+S[104])
ogrxn[74] = kcat[74]*E[74] * S[91]/(Km[91,74]+S[91]) * S[106]/(Km[106,74]+S[106])
ogrxn[75] = kcat[75]*E[75] * S[16]/(Km[16,75]+S[16])
ogrxn[76] = kcat[76]*E[76]
ogrxn[77] = kcat[77]*E[77] * S[18]/(Km[18,77]+S[18])
ogrxn[78] = kcat[78]*E[78] * S[15]/(Km[15,78]+S[15]) * S[31]/(Km[31,78]+S[31])
ogrxn[79] = kcat[79]*E[79] * S[15]/(Km[15,79]+S[15]) * S[31]/(Km[31,79]+S[31])
ogrxn[80] = kcat[80]*E[80] * S[16]/(Km[16,80]+S[16]) * S[105]/(Km[105,80]+S[105])
ogrxn[81] = kcat[81]*E[81] * S[104]/(Km[104,81]+S[104])
ogrxn[82] = kcat[82]*E[82] * S[105]/(Km[105,82]+S[105]) * S[122]/(Km[122,82]+S[122])
ogrxn[83] = kcat[83]*E[83] * S[86]/(Km[86,83]+S[86]) * S[104]/(Km[104,83]+S[104])
ogrxn[84] = kcat[84]*E[84] * S[122]/(Km[122,84]+S[122])
ogrxn[85] = kcat[85]*E[85] * S[64]/(Km[64,85]+S[64]) * S[122]/(Km[122,85]+S[122])
ogrxn[86] = kcat[86]*E[86] * S[21]/(Km[21,86]+S[21]) * S[22]/(Km[22,86]+S[22])
ogrxn[87] = kcat[87]*E[87] * S[16]/(Km[16,87]+S[16]) * S[29]/(Km[29,87]+S[29]) * S[31]/(Km[31,87]+S[31]) * S[64]/(Km[64,87]+S[64]) * S[107]/(Km[107,87]+S[107])
ogrxn[88] = kcat[88]*E[88] * S[64]/(Km[64,88]+S[64]) * S[110]/(Km[110,88]+S[110])
ogrxn[89] = kcat[89]*E[89] * S[29]/(Km[29,89]+S[29]) * S[31]/(Km[31,89]+S[31]) * S[62]/(Km[62,89]+S[62])
ogrxn[90] = kcat[90]*E[90] * S[29]/(Km[29,90]+S[29]) * S[31]/(Km[31,90]+S[31])
ogrxn[91] = kcat[91]*E[91] * S[16]/(Km[16,91]+S[16]) * S[129]/(Km[129,91]+S[129])
ogrxn[92] = kcat[92]*E[92] * S[21]/(Km[21,92]+S[21]) * S[62]/(Km[62,92]+S[62]) * S[107]/(Km[107,92]+S[107])
ogrxn[93] = kcat[93]*E[93] * S[21]/(Km[21,93]+S[21]) * S[107]/(Km[107,93]+S[107])
ogrxn[94] = kcat[94]*E[94] * S[64]/(Km[64,94]+S[64]) * S[106]/(Km[106,94]+S[106])
ogrxn[95] = kcat[95]*E[95] * S[31]/(Km[31,95]+S[31]) * S[64]/(Km[64,95]+S[64])
ogrxn[96] = kcat[96]*E[96] * S[129]/(Km[129,96]+S[129])
ogrxn[97] = kcat[97]*E[97] * S[31]/(Km[31,97]+S[31]) * S[62]/(Km[62,97]+S[62]) * S[104]/(Km[104,97]+S[104]) * S[125]/(Km[125,97]+S[125])
ogrxn[98] = kcat[98]*E[98] * S[64]/(Km[64,98]+S[64]) * S[107]/(Km[107,98]+S[107]) * S[122]/(Km[122,98]+S[122]) * S[134]/(Km[134,98]+S[134])
ogrxn[99] = kcat[99]*E[99] * S[16]/(Km[16,99]+S[16]) * S[64]/(Km[64,99]+S[64]) * S[104]/(Km[104,99]+S[104]) * S[107]/(Km[107,99]+S[107]) * S[122]/(Km[122,99]+S[122])
ogrxn[100] = kcat[100]*E[100] * S[29]/(Km[29,100]+S[29]) * S[31]/(Km[31,100]+S[31]) * S[64]/(Km[64,100]+S[64]) * S[107]/(Km[107,100]+S[107]) * S[122]/(Km[122,100]+S[122])
ogrxn[101] = kcat[101]*E[101] * S[29]/(Km[29,101]+S[29]) * S[31]/(Km[31,101]+S[31]) * S[42]/(Km[42,101]+S[42]) * S[107]/(Km[107,101]+S[107])
ogrxn[102] = kcat[102]*E[102] * S[64]/(Km[64,102]+S[64])
ogrxn[103] = kcat[103]*E[103] * S[31]/(Km[31,103]+S[31]) * S[64]/(Km[64,103]+S[64]) * S[107]/(Km[107,103]+S[107])
ogrxn[104] = kcat[104]*E[104] * S[6]/(Km[6,104]+S[6]) * S[64]/(Km[64,104]+S[64]) * S[104]/(Km[104,104]+S[104])
ogrxn[105] = kcat[105]*E[105] * S[29]/(Km[29,105]+S[29]) * S[31]/(Km[31,105]+S[31]) * S[107]/(Km[107,105]+S[107])
ogrxn[106] = kcat[106]*E[106] * S[31]/(Km[31,106]+S[31]) * S[62]/(Km[62,106]+S[62]) * S[125]/(Km[125,106]+S[125]) * S[129]/(Km[129,106]+S[129])
ogrxn[107] = kcat[107]*E[107] * S[64]/(Km[64,107]+S[64]) * S[104]/(Km[104,107]+S[104])
ogrxn[108] = kcat[108]*E[108] * S[64]/(Km[64,108]+S[64]) * S[107]/(Km[107,108]+S[107]) * S[122]/(Km[122,108]+S[122])
ogrxn[109] = kcat[109]*E[109] * S[21]/(Km[21,109]+S[21]) * S[25]/(Km[25,109]+S[25]) * S[104]/(Km[104,109]+S[104])
ogrxn[110] = kcat[110]*E[110] * S[29]/(Km[29,110]+S[29])
ogrxn[111] = kcat[111]*E[111] * S[24]/(Km[24,111]+S[24]) * S[27]/(Km[27,111]+S[27])
ogrxn[112] = kcat[112]*E[112] * S[16]/(Km[16,112]+S[16]) * S[68]/(Km[68,112]+S[68])
ogrxn[113] = kcat[113]*E[113] * S[104]/(Km[104,113]+S[104])
ogrxn[114] = kcat[114]*E[114] * S[129]/(Km[129,114]+S[129])
ogrxn[115] = kcat[115]*E[115] * S[104]/(Km[104,115]+S[104]) * S[118]/(Km[118,115]+S[118])
ogrxn[116] = kcat[116]*E[116] * S[104]/(Km[104,116]+S[104]) * S[134]/(Km[134,116]+S[134])
ogrxn[117] = kcat[117]*E[117] * S[104]/(Km[104,117]+S[104]) * S[134]/(Km[134,117]+S[134])
ogrxn[118] = kcat[118]*E[118] * S[18]/(Km[18,118]+S[18]) * S[134]/(Km[134,118]+S[134])
ogrxn[119] = kcat[119]*E[119] * S[136]/(Km[136,119]+S[136])
ogrxn[120] = kcat[120]*E[120] * S[42]/(Km[42,120]+S[42])
ogrxn[121] = kcat[121]*E[121] * S[89]/(Km[89,121]+S[89])
ogrxn[122] = kcat[122]*E[122] * S[62]/(Km[62,122]+S[62])
ogrxn[123] = kcat[123]*E[123] * S[64]/(Km[64,123]+S[64])
ogrxn[124] = kcat[124]*E[124] * S[21]/(Km[21,124]+S[21]) * S[104]/(Km[104,124]+S[104])
ogrxn[125] = kcat[125]*E[125] * S[21]/(Km[21,125]+S[21]) * S[106]/(Km[106,125]+S[106])
ogrxn[126] = kcat[126]*E[126] * S[31]/(Km[31,126]+S[31]) * S[46]/(Km[46,126]+S[46]) * S[107]/(Km[107,126]+S[107]) * S[113]/(Km[113,126]+S[113])
ogrxn[127] = kcat[127]*E[127] * S[72]/(Km[72,127]+S[72])
ogrxn[128] = kcat[128]*E[128] * S[62]/(Km[62,128]+S[62])
ogrxn[129] = kcat[129]*E[129]
ogrxn[130] = kcat[130]*E[130]
ogrxn[131] = kcat[131]*E[131] * S[31]/(Km[31,131]+S[31]) * S[64]/(Km[64,131]+S[64])
ogrxn[132] = kcat[132]*E[132] * S[107]/(Km[107,132]+S[107])
ogrxn[133] = kcat[133]*E[133] * S[68]/(Km[68,133]+S[68]) * S[104]/(Km[104,133]+S[104])
ogrxn[134] = kcat[134]*E[134] * S[105]/(Km[105,134]+S[105])
ogrxn[135] = kcat[135]*E[135] * S[106]/(Km[106,135]+S[106])
ogrxn[136] = kcat[136]*E[136] * S[107]/(Km[107,136]+S[107])
ogrxn[137] = kcat[137]*E[137]
ogrxn[138] = kcat[138]*E[138]
ogrxn[139] = kcat[139]*E[139] * S[105]/(Km[105,139]+S[105])
ogrxn[140] = kcat[140]*E[140] * S[107]/(Km[107,140]+S[107])
ogrxn[141] = kcat[141]*E[141] * S[31]/(Km[31,141]+S[31]) * S[125]/(Km[125,141]+S[125])
ogrxn[142] = kcat[142]*E[142] * S[31]/(Km[31,142]+S[31]) * S[62]/(Km[62,142]+S[62])
ogrxn[143] = kcat[143]*E[143] * S[29]/(Km[29,143]+S[29])
ogrxn[144] = kcat[144]*E[144]
ogrxn[145] = kcat[145]*E[145]
ogrxn[146] = kcat[146]*E[146] * S[31]/(Km[31,146]+S[31]) * S[142]/(Km[142,146]+S[142])
ogrxn[147] = kcat[147]*E[147] * S[31]/(Km[31,147]+S[31]) * S[62]/(Km[62,147]+S[62]) * S[142]/(Km[142,147]+S[142])
ogrxn[148] = kcat[148]*E[148] * S[62]/(Km[62,148]+S[62])
ogrxn[149] = kcat[149]*E[149] * S[31]/(Km[31,149]+S[31]) * S[68]/(Km[68,149]+S[68])
ogrxn[150] = kcat[150]*E[150]
ogrxn[151] = kcat[151]*E[151] * S[31]/(Km[31,151]+S[31]) * S[62]/(Km[62,151]+S[62])
ogrxn[152] = kcat[152]*E[152] * S[31]/(Km[31,152]+S[31])
ogrxn[153] = kcat[153]*E[153] * S[31]/(Km[31,153]+S[31])
ogrxn[154] = kcat[154]*E[154] * S[29]/(Km[29,154]+S[29]) * S[31]/(Km[31,154]+S[31])
ogrxn[155] = kcat[155]*E[155]
ogrxn[156] = kcat[156]*E[156]
ogrxn[157] = kcat[157]*E[157]
ogrxn[158] = kcat[158]*E[158] * S[29]/(Km[29,158]+S[29]) * S[72]/(Km[72,158]+S[72])
ogrxn[159] = kcat[159]*E[159] * S[104]/(Km[104,159]+S[104])
ogrxn[160] = kcat[160]*E[160] * S[31]/(Km[31,160]+S[31]) * S[62]/(Km[62,160]+S[62])
ogrxn[161] = kcat[161]*E[161] * S[31]/(Km[31,161]+S[31])
ogrxn[162] = kcat[162]*E[162] * S[142]/(Km[142,162]+S[142])
ogrxn[163] = kcat[163]*E[163] * S[41]/(Km[41,163]+S[41])
ogrxn[164] = kcat[164]*E[164] * S[72]/(Km[72,164]+S[72])
ogrxn[165] = kcat[165]*E[165] * S[31]/(Km[31,165]+S[31])
ogrxn[166] = kcat[166]*E[166] * S[142]/(Km[142,166]+S[142])
ogrxn[167] = kcat[167]*E[167] * S[41]/(Km[41,167]+S[41])
ogrxn[168] = kcat[168]*E[168] * S[72]/(Km[72,168]+S[72])
ogrxn[169] = kcat[169]*E[169] * S[31]/(Km[31,169]+S[31]) * S[140]/(Km[140,169]+S[140])
ogrxn[170] = kcat[170]*E[170] * S[31]/(Km[31,170]+S[31]) * S[34]/(Km[34,170]+S[34])
ogrxn[171] = kcat[171]*E[171] * S[31]/(Km[31,171]+S[31]) * S[60]/(Km[60,171]+S[60])
ogrxn[172] = kcat[172]*E[172] * S[31]/(Km[31,172]+S[31]) * S[141]/(Km[141,172]+S[141])
ogrxn[173] = kcat[173]*E[173] * S[31]/(Km[31,173]+S[31]) * S[38]/(Km[38,173]+S[38])
ogrxn[174] = kcat[174]*E[174] * S[31]/(Km[31,174]+S[31]) * S[71]/(Km[71,174]+S[71])
ogrxn[175] = kcat[175]*E[175] * S[24]/(Km[24,175]+S[24]) * S[31]/(Km[31,175]+S[31])
RX_conc = TXTL[:RX_concentration]
RL_conc = TXTL[:RL_concentration]
  
ogrxn[176] = VX/RX_conc*length_factor_transcription*(S[1]/(KX*X_tau_factor+(1+X_tau_factor)*S[1])) #mM/hr
ogrxn[177] = VX/RX_conc*length_factor_transcription*(S[1]/(KX*X_tau_factor+(1+X_tau_factor)*S[1])) * S[31]/(km_tx+S[31]) * S[41]/(km_tx+S[41]) * S[72]/(km_tx+S[72]) * S[142]/(km_tx+S[142])
ogrxn[178] = S[152] * mRNA_degradation_rate
ogrxn[179] = VL/RL_conc*length_factor_translation*(S[152]/(KL*L_tau_factor+(1+L_tau_factor)*S[152]))
ogrxn[180] = VL/RL_conc*length_factor_translation*(S[152]/(KL*L_tau_factor+(1+L_tau_factor)*S[152])) * S[72]/(0.05+S[72])

#= ogrxn[176] = kcat[176]*E[176] * S[1]/(Km[1,176]+S[1]) * S[151]/(Km[151,176]+S[151])
ogrxn[177] = kcat[177]*E[177] * S[31]/(Km[31,177]+S[31]) * S[41]/(Km[41,177]+S[41]) * S[72]/(Km[72,177]+S[72]) * S[74]/(Km[74,177]+S[74]) * S[142]/(Km[142,177]+S[142]) * S[147]/(Km[147,177]+S[147])
ogrxn[178] = kcat[178]*E[178] * S[152]/(Km[152,178]+S[152])
ogrxn[179] = kcat[179]*E[179] * S[149]/(Km[149,179]+S[149]) * S[152]/(Km[152,179]+S[152])
ogrxn[180] = kcat[180]*E[180] * S[23]/(Km[23,180]+S[23]) * S[26]/(Km[26,180]+S[26]) * S[28]/(Km[28,180]+S[28]) * S[30]/(Km[30,180]+S[30]) * S[43]/(Km[43,180]+S[43]) * S[63]/(Km[63,180]+S[63]) * S[65]/(Km[65,180]+S[65]) * S[69]/(Km[69,180]+S[69]) * S[72]/(Km[72,180]+S[72]) * S[74]/(Km[74,180]+S[74]) * S[79]/(Km[79,180]+S[79]) * S[83]/(Km[83,180]+S[83]) * S[88]/(Km[88,180]+S[88]) * S[90]/(Km[90,180]+S[90]) * S[98]/(Km[98,180]+S[98]) * S[115]/(Km[115,180]+S[115]) * S[119]/(Km[119,180]+S[119]) * S[130]/(Km[130,180]+S[130]) * S[135]/(Km[135,180]+S[135]) * S[137]/(Km[137,180]+S[137]) * S[139]/(Km[139,180]+S[139]) * S[144]/(Km[144,180]+S[144]) * S[150]/(Km[150,180]+S[150])
 =#
 ogrxn[181] = kcat[181]*E[181] * S[22]/(Km[22,181]+S[22]) * S[31]/(Km[31,181]+S[31])
 ogrxn[182] = kcat[182]*E[182] * S[25]/(Km[25,182]+S[25]) * S[31]/(Km[31,182]+S[31])
 ogrxn[183] = kcat[183]*E[183] * S[27]/(Km[27,183]+S[27]) * S[31]/(Km[31,183]+S[31])
 ogrxn[184] = kcat[184]*E[184] * S[29]/(Km[29,184]+S[29]) * S[31]/(Km[31,184]+S[31]) 
 ogrxn[185] = kcat[185]*E[185] * S[31]/(Km[31,185]+S[31]) * S[42]/(Km[42,185]+S[42]) 
 ogrxn[186] = kcat[186]*E[186] * S[31]/(Km[31,186]+S[31]) * S[64]/(Km[64,186]+S[64])
 ogrxn[187] = kcat[187]*E[187] * S[31]/(Km[31,187]+S[31]) * S[62]/(Km[62,187]+S[62])
 ogrxn[188] = kcat[188]*E[188] * S[31]/(Km[31,188]+S[31]) * S[68]/(Km[68,188]+S[68]) 
 ogrxn[189] = kcat[189]*E[189] * S[31]/(Km[31,189]+S[31]) * S[78]/(Km[78,189]+S[78]) 
 ogrxn[190] = kcat[190]*E[190] * S[31]/(Km[31,190]+S[31]) * S[82]/(Km[82,190]+S[82]) 
 ogrxn[191] = kcat[191]*E[191] * S[31]/(Km[31,191]+S[31]) * S[87]/(Km[87,191]+S[87]) 
 ogrxn[192] = kcat[192]*E[192] * S[31]/(Km[31,192]+S[31]) * S[89]/(Km[89,192]+S[89]) 
 ogrxn[193] = kcat[193]*E[193] * S[31]/(Km[31,193]+S[31]) * S[97]/(Km[97,193]+S[97]) 
 ogrxn[194] = kcat[194]*E[194] * S[31]/(Km[31,194]+S[31]) * S[114]/(Km[114,194]+S[114]) 
 ogrxn[195] = kcat[195]*E[195] * S[31]/(Km[31,195]+S[31]) * S[118]/(Km[118,195]+S[118]) 
 ogrxn[196] = kcat[196]*E[196] * S[31]/(Km[31,196]+S[31]) * S[129]/(Km[129,196]+S[129]) 
 ogrxn[197] = kcat[197]*E[197] * S[31]/(Km[31,197]+S[31]) * S[134]/(Km[134,197]+S[134]) 
 ogrxn[198] = kcat[198]*E[198] * S[31]/(Km[31,198]+S[31]) * S[136]/(Km[136,198]+S[136])
 ogrxn[199] = kcat[199]*E[199] * S[31]/(Km[31,199]+S[31]) * S[138]/(Km[138,199]+S[138]) 
 ogrxn[200] = kcat[200]*E[200] * S[31]/(Km[31,200]+S[31]) * S[143]/(Km[143,200]+S[143])

return ogrxn

end

#= using ForwardDiff
S = ones(157);
jac = ForwardDiff.jacobian(ogrxnsetup!,S);
jac = jac';
nz = findall(x->x!=0&&!isnan(x),jac);
jacnz = jac[nz];
 =#
# writedlm("jacobian_gammas.txt",jac)

