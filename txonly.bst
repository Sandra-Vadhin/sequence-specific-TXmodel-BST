// s-law model for degfp transcription in mytxtl -

#static::start
E_gp,E_pgm,E_glk_atp,E_pgi,E_pfk,E_fdp,E_fbaA,E_tpiA,E_gpsA,E_gapA,E_pgk,E_gpm,E_eno,E_pyk,E_pck,E_ppc,E_pdh,E_pps,E_zwf,E_pgl,E_gnd,E_rpe,E_rpi,E_talAB,E_tkt,E_edd,E_eda,E_gltA,E_acn,E_icd,E_sucAB,E_sucCD,E_sdh,E_frd,E_fum,E_mdh,E_cyco,E_atp,E_DNP_leak,E_nuo,E_pnt1,E_pnt2,E_ppa,E_aceA,E_aceB,E_maeA,E_maeB,E_pta,E_ackA,E_acs,E_adhE,E_ldh,E_pflAB,E_alaAC,E_arg,E_aspC,E_asnB,E_asnA,E_cysEMK,E_gltBD,E_gdhA,E_glnA,E_glyA,E_his,E_ile,E_leu,E_lys,E_met,E_phe,E_pro,E_serABC,E_thr,E_trp,E_tyr,E_val,E_arg_deg,E_asp_deg,E_gly_deg,E_mglx_deg,E_ser_deg,E_pro_deg,E_thr_deg1,E_thr_deg2,E_thr_deg3,E_trp_deg,E_cys_deg,E_lys_deg,E_gln_deg,E_glu_deg,E_gaba_deg,E_chor,E_fol_e,E_fol_1,E_fol_2a,E_fol_2b,E_fol_3,E_fol_4,E_gly_fol,E_mthfd,E_mthfc,E_mthfr2,E_prpp_syn,E_or_syn_1,E_or_syn_2,E_omp_syn,E_ump_syn,E_ctp,E_A_syn_1,E_A_syn_2,E_A_syn_3,E_A_syn_4,E_A_syn_5,E_A_syn_6,E_A_syn_7,E_A_syn_8,E_A_syn_9,E_A_syn_10,E_A_syn_12,E_xmp_syn,E_gmp_syn,E_NDPase,E_atp_adp,E_ADPase,E_ndk,E_atp_ump,E_atp_cmp,E_atp_gmp,E_adk_atp,E_tRNA_c_charging_M_ala_L_c,E_tRNA_c_charging_M_arg_L_c,E_tRNA_c_charging_M_asn_L_c,E_tRNA_c_charging_M_asp_L_c,E_tRNA_c_charging_M_cys_L_c,E_tRNA_c_charging_M_glu_L_c,E_tRNA_c_charging_M_gln_L_c,E_tRNA_c_charging_M_gly_L_c,E_tRNA_c_charging_M_his_L_c,E_tRNA_c_charging_M_ile_L_c,E_tRNA_c_charging_M_leu_L_c,E_tRNA_c_charging_M_lys_L_c,E_tRNA_c_charging_M_met_L_c,E_tRNA_c_charging_M_phe_L_c,E_tRNA_c_charging_M_pro_L_c,E_tRNA_c_charging_M_ser_L_c,E_tRNA_c_charging_M_thr_L_c,E_tRNA_c_charging_M_trp_L_c,E_tRNA_c_charging_M_tyr_L_c,E_tRNA_c_charging_M_val_L_c
#static::end

#dynamic::start
GENE_gfp,M_10fthf_c,M_13dpg_c,M_2ddg6p_c,M_2pg_c,M_3pg_c,M_4abz_c,M_4adochor_c,M_5mthf_c,M_5pbdra,M_6pgc_c,M_6pgl_c,M_78dhf_c,M_78mdp_c,M_ac_c,M_accoa_c,M_actp_c,M_adp_c,M_aicar_c,M_air_c,M_akg_c,M_ala_L_c,M_ala_L_c_tRNA_c,M_amp_c,M_arg_L_c,M_arg_L_c_tRNA_c,M_asn_L_c,M_asn_L_c_tRNA_c,M_asp_L_c,M_asp_L_c_tRNA_c,M_atp_c,M_cadav_c,M_cair_c,M_cdp_c,M_chor_c,M_cit_c,M_clasp_c,M_cmp_c,M_co2_c,M_coa_c,M_ctp_c,M_cys_L_c,M_cys_L_c_tRNA_c,M_dhap_c,M_dhf_c,M_e4p_c,M_etoh_c,M_f6p_c,M_faicar_c,M_fdp_c,M_fgam_c,M_fgar_c,M_for_c,M_fum_c,M_g1p_c,M_g3p_c,M_g6p_c,M_gaba_c,M_gar_c,M_gdp_c,M_glc_D_c,M_gln_L_c,M_gln_L_c_tRNA_c,M_glu_L_c,M_glu_L_c_tRNA_c,M_glx_c,M_gly3p_c,M_gly_L_c,M_gly_L_c_tRNA_c,M_glycoA_c,M_gmp_c,M_gtp_c,M_h2o2_c,M_h2o_c,M_h2s_c,M_h_c,M_hco3_c,M_his_L_c,M_his_L_c_tRNA_c,M_hx_c,M_icit_c,M_ile_L_c,M_ile_L_c_tRNA_c,M_imp_c,M_indole_c,M_lac_D_c,M_leu_L_c,M_leu_L_c_tRNA_c,M_lys_L_c,M_lys_L_c_tRNA_c,M_mal_L_c,M_maltodextrin3_c,M_maltodextrin4_c,M_maltodextrin5_c,M_maltodextrin6_c,M_maltose_c,M_met_L_c,M_met_L_c_tRNA_c,M_methf_c,M_mglx_c,M_mlthf_c,M_mql8_c,M_mqn8_c,M_nad_c,M_nadh_c,M_nadp_c,M_nadph_c,M_nh3_c,M_o2_c,M_oaa_c,M_omp_c,M_or_c,M_pep_c,M_phe_L_c,M_phe_L_c_tRNA_c,M_pi_c,M_ppi_c,M_pro_L_c,M_pro_L_c_tRNA_c,M_prop_c,M_prpp_c,M_pyr_c,M_q8_c,M_q8h2_c,M_r5p_c,M_ru5p_D_c,M_s7p_c,M_saicar_c,M_ser_L_c,M_ser_L_c_tRNA_c,M_succ_c,M_succoa_c,M_thf_c,M_thr_L_c,M_thr_L_c_tRNA_c,M_trp_L_c,M_trp_L_c_tRNA_c,M_tyr_L_c,M_tyr_L_c_tRNA_c,M_udp_c,M_ump_c,M_utp_c,M_val_L_c,M_val_L_c_tRNA_c,M_xmp_c,M_xu5p_D_c,OPEN_GENE_gfp,PROTEIN_gfp,RIBOSOME,RIBOSOME_START_gfp,RNAP,mRNA_gfp,tRNA_c,M_co2_e,M_h_e,M_o2_e,gfp_e,INDUCER
#dynamic::end

#structure::start
    // name::{reactants} --> {products}
    // 1. name: reaction name (used in rate section)
    // 2. {reactants}: list of dynamic species acting as a reactants
    // 3. {products}: list of dynamic species acting as a reactants
transcriptional_initiation_gfp::{GENE_gfp,RNAP,INDUCER} --> {OPEN_GENE_gfp}
transcription_gfp::{OPEN_GENE_gfp,M_gtp_c,M_ctp_c,M_utp_c,M_atp_c,M_h2o_c} --> {mRNA_gfp,GENE_gfp,RNAP,M_ppi_c,INDUCER}
mRNA_degradation_gfp::{mRNA_gfp} --> {M_gmp_c,M_cmp_c,M_ump_c,M_amp_c}
R_ppa::{M_ppi_c,M_h2o_c} --> {M_pi_c}
//tRNA_c_c_exchange::{} --> {tRNA_c}
//tRNA_c_c_exchange_reverse::{tRNA_c} --> {}
M_h2o_c_exchange::{} --> {M_h2o_c}
M_h2o_c_exchange_reverse::{M_h2o_c} --> {}
M_pi_c_exchange::{} --> {M_pi_c}
M_pi_c_exchange_reverse::{M_pi_c} --> {}
Exchange_M_atp_c::{} --> {M_atp_c}
Exchange_M_ctp_c::{} --> {M_ctp_c}
Exchange_M_gtp_c::{} --> {M_gtp_c}
Exchange_M_utp_c::{} --> {M_utp_c}
Exchange_M_atp_c_REV::{M_atp_c} --> {}
Exchange_M_ctp_c_REV::{M_ctp_c} --> {}
Exchange_M_gtp_c_REV::{M_gtp_c} --> {}
Exchange_M_utp_c_REV::{M_utp_c} --> {}
Exchange_M_amp_c::{M_amp_c} --> {}
Exchange_M_cmp_c::{M_cmp_c} --> {}
Exchange_M_gmp_c::{M_gmp_c} --> {}
Exchange_M_ump_c::{M_ump_c} --> {}
//Exchange_M_gdp_c::{M_gdp_c} --> {}
Exchange_M_amp_c_REV::{} --> {M_amp_c}
Exchange_M_cmp_c_REV::{} --> {M_cmp_c}
Exchange_M_gmp_c_REV::{} --> {M_gmp_c}
Exchange_M_ump_c_REV::{} --> {M_ump_c}
//Exchange_M_gdp_c_REV::{} --> {M_gdp_c}
#structure::end

#rate::start
transcriptional_initiation_gfp::{GENE_gfp,RNAP,INDUCER}
transcription_gfp::{OPEN_GENE_gfp,M_gtp_c,M_ctp_c,M_utp_c,M_atp_c,M_h2o_c}
mRNA_degradation_gfp::{mRNA_gfp}
//translation_initiation_gfp::{mRNA_gfp,RIBOSOME}
//translation_gfp::{RIBOSOME_START_gfp,M_gtp_c,M_h2o_c,M_aa_tRNA_c}
R_ppa::{E_ppa,M_ppi_c,M_h2o_c}
//tRNA_c_c_exchange::{}
//tRNA_c_c_exchange_reverse::{tRNA_c}
M_h2o_c_exchange::{}
M_h2o_c_exchange_reverse::{M_h2o_c}
M_pi_c_exchange::{}
M_pi_c_exchange_reverse::{M_pi_c}
Exchange_M_atp_c::{}
Exchange_M_ctp_c::{}
Exchange_M_gtp_c::{}
Exchange_M_utp_c::{}
Exchange_M_amp_c::{M_amp_c}
Exchange_M_cmp_c::{M_cmp_c}
Exchange_M_gmp_c::{M_gmp_c}
Exchange_M_ump_c::{M_ump_c}
//Exchange_M_gdp_c::{M_gdp_c}
Exchange_M_atp_c_REV::{M_atp_c}
Exchange_M_ctp_c_REV::{M_ctp_c}
Exchange_M_gtp_c_REV::{M_gtp_c}
Exchange_M_utp_c_REV::{M_utp_c}
Exchange_M_amp_c_REV::{}
Exchange_M_cmp_c_REV::{}
Exchange_M_gmp_c_REV::{}
Exchange_M_ump_c_REV::{}
//Exchange_M_gdp_c_REV::{}
//Exchange_aa_c::{}
//Exchange_aa_c_REV::{M_aa_c}
#rate::end

#stoichiometry::start
    // rxn_name,species,coeff
    // 1. rxn_name: reaction name (used in rate section)
    // 2. species: species in reaction
    // 3. coeff: stoichiometric coefficient value of species in reaction
transcription_gfp,M_gtp_c,195.0
transcription_gfp,M_ctp_c,157.0
transcription_gfp,M_utp_c,157.0
transcription_gfp,M_atp_c,208.0
transcription_gfp,M_h2o_c,717.0
transcription_gfp,M_ppi_c,717.0
mRNA_degradation_gfp,M_gmp_c,195.0
mRNA_degradation_gfp,M_cmp_c,157.0
mRNA_degradation_gfp,M_ump_c,157.0
mRNA_degradation_gfp,M_amp_c,208.0
//translation_gfp,M_gtp_c,478.0
//translation_gfp,M_h2o_c,478.0
//translation_gfp,M_gdp_c,478.0
//translation_gfp,M_pi_c,478.0
//translation_gfp,tRNA_c,239.0
//translation_gfp,M_aa_c,239.0
R_ppa,M_pi_c,2.0
#stoichiometry::end