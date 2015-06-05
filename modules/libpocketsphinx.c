/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */


/**
 * @file acmod.c Acoustic model structures for PocketSphinx.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

/* System headers. */
#include <assert.h>
#include <string.h>

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>
#include <sphinxbase/err.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/byteorder.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/bio.h>

/* Local headers. */
#include "cmdln_macro.h"
#include "acmod.h"
#include "s2_semi_mgau.h"
#include "ptm_mgau.h"
#include "ms_mgau.h"

/* Feature and front-end parameters that may be in feat.params */
static const arg_t feat_defn[] = {
    waveform_to_cepstral_command_line_macro(),
    cepstral_to_feature_command_line_macro(),
    CMDLN_EMPTY_OPTION
};

#ifndef WORDS_BIGENDIAN
#define WORDS_BIGENDIAN 1
#endif

static int32 acmod_process_mfcbuf(acmod_t *acmod);

static int
acmod_init_am(acmod_t *acmod)
{
    char const *mdeffn, *tmatfn, *mllrfn, *hmmdir;

    /* Read model definition. */
    if ((mdeffn = cmd_ln_str_r(acmod->config, "-mdef")) == NULL) {
	if ((hmmdir = cmd_ln_str_r(acmod->config, "-hmm")) == NULL) {
	    E_ERROR("Acoustic model definition is not specified neither with -mdef option nor with -hmm\n");
	} else {
	    E_ERROR("Folder '%s' does not contain acoustic model definition 'mdef'\n", hmmdir);
	}
        return -1;
    }

    if ((acmod->mdef = bin_mdef_read(acmod->config, mdeffn)) == NULL) {
        E_ERROR("Failed to read acoustic model definition from %s\n", mdeffn);
        return -1;
    }

    /* Read transition matrices. */
    if ((tmatfn = cmd_ln_str_r(acmod->config, "-tmat")) == NULL) {
        E_ERROR("No tmat file specified\n");
        return -1;
    }
    acmod->tmat = tmat_init(tmatfn, acmod->lmath,
                            cmd_ln_float32_r(acmod->config, "-tmatfloor"),
                            TRUE);

    /* Read the acoustic models. */
    if ((cmd_ln_str_r(acmod->config, "-mean") == NULL)
        || (cmd_ln_str_r(acmod->config, "-var") == NULL)
        || (cmd_ln_str_r(acmod->config, "-tmat") == NULL)) {
        E_ERROR("No mean/var/tmat files specified\n");
        return -1;
    }

    if (cmd_ln_str_r(acmod->config, "-senmgau")) {
        E_INFO("Using general multi-stream GMM computation\n");
        acmod->mgau = ms_mgau_init(acmod, acmod->lmath, acmod->mdef);
        if (acmod->mgau == NULL)
            return -1;
    }
    else {
        E_INFO("Attempting to use SCHMM computation module\n");
        if ((acmod->mgau = s2_semi_mgau_init(acmod)) == NULL) {
            E_INFO("Attempting to use PTHMM computation module\n");
            if ((acmod->mgau = ptm_mgau_init(acmod, acmod->mdef)) == NULL) {
                E_INFO("Falling back to general multi-stream GMM computation\n");
                acmod->mgau = ms_mgau_init(acmod, acmod->lmath, acmod->mdef);
                if (acmod->mgau == NULL)
                    return -1;
            }
        }
    }

    /* If there is an MLLR transform, apply it. */
    if ((mllrfn = cmd_ln_str_r(acmod->config, "-mllr"))) {
        ps_mllr_t *mllr = ps_mllr_read(mllrfn);
        if (mllr == NULL)
            return -1;
        acmod_update_mllr(acmod, mllr);
    }

    return 0;
}

static int
acmod_init_feat(acmod_t *acmod)
{
    acmod->fcb = 
        feat_init(cmd_ln_str_r(acmod->config, "-feat"),
                  cmn_type_from_str(cmd_ln_str_r(acmod->config,"-cmn")),
                  cmd_ln_boolean_r(acmod->config, "-varnorm"),
                  agc_type_from_str(cmd_ln_str_r(acmod->config, "-agc")),
                  1, cmd_ln_int32_r(acmod->config, "-ceplen"));
    if (acmod->fcb == NULL)
        return -1;

    if (cmd_ln_str_r(acmod->config, "-lda")) {
        E_INFO("Reading linear feature transformation from %s\n",
               cmd_ln_str_r(acmod->config, "-lda"));
        if (feat_read_lda(acmod->fcb,
                          cmd_ln_str_r(acmod->config, "-lda"),
                          cmd_ln_int32_r(acmod->config, "-ldadim")) < 0)
            return -1;
    }

    if (cmd_ln_str_r(acmod->config, "-svspec")) {
        int32 **subvecs;
        E_INFO("Using subvector specification %s\n", 
               cmd_ln_str_r(acmod->config, "-svspec"));
        if ((subvecs = parse_subvecs(cmd_ln_str_r(acmod->config, "-svspec"))) == NULL)
            return -1;
        if ((feat_set_subvecs(acmod->fcb, subvecs)) < 0)
            return -1;
    }

    if (cmd_ln_exists_r(acmod->config, "-agcthresh")
        && 0 != strcmp(cmd_ln_str_r(acmod->config, "-agc"), "none")) {
        agc_set_threshold(acmod->fcb->agc_struct,
                          cmd_ln_float32_r(acmod->config, "-agcthresh"));
    }

    if (acmod->fcb->cmn_struct
        && cmd_ln_exists_r(acmod->config, "-cmninit")) {
        char *c, *cc, *vallist;
        int32 nvals;

        vallist = ckd_salloc(cmd_ln_str_r(acmod->config, "-cmninit"));
        c = vallist;
        nvals = 0;
        while (nvals < acmod->fcb->cmn_struct->veclen
               && (cc = strchr(c, ',')) != NULL) {
            *cc = '\0';
            acmod->fcb->cmn_struct->cmn_mean[nvals] = FLOAT2MFCC(atof(c));
            c = cc + 1;
            ++nvals;
        }
        if (nvals < acmod->fcb->cmn_struct->veclen && *c != '\0') {
            acmod->fcb->cmn_struct->cmn_mean[nvals] = FLOAT2MFCC(atof(c));
        }
        ckd_free(vallist);
    }
    return 0;
}

int
acmod_fe_mismatch(acmod_t *acmod, fe_t *fe)
{
    /* Output vector dimension needs to be the same. */
    if (cmd_ln_int32_r(acmod->config, "-ceplen") != fe_get_output_size(fe)) {
	E_ERROR("Configured feature length %d doesn't match feature extraction output size %d\n", 
		cmd_ln_int32_r(acmod->config, "-ceplen"), 
		fe_get_output_size(fe));
        return TRUE;
    }
    /* Feature parameters need to be the same. */
    /* ... */
    return FALSE;
}

int
acmod_feat_mismatch(acmod_t *acmod, feat_t *fcb)
{
    /* Feature type needs to be the same. */
    if (0 != strcmp(cmd_ln_str_r(acmod->config, "-feat"), feat_name(fcb)))
        return TRUE;
    /* Input vector dimension needs to be the same. */
    if (cmd_ln_int32_r(acmod->config, "-ceplen") != feat_cepsize(fcb))
        return TRUE;
    /* FIXME: Need to check LDA and stuff too. */
    return FALSE;
}

acmod_t *
acmod_init(cmd_ln_t *config, logmath_t *lmath, fe_t *fe, feat_t *fcb)
{
    acmod_t *acmod;
    char const *featparams;

    acmod = ckd_calloc(1, sizeof(*acmod));
    acmod->config = cmd_ln_retain(config);
    acmod->lmath = lmath;
    acmod->state = ACMOD_IDLE;

    /* Look for feat.params in acoustic model dir. */
    if ((featparams = cmd_ln_str_r(acmod->config, "-featparams"))) {
        if (cmd_ln_parse_file_r(acmod->config, feat_defn, featparams, FALSE) != NULL) {
	    E_INFO("Parsed model-specific feature parameters from %s\n", featparams);
        }
    }

    /* Initialize feature computation. */
    if (fe) {
        if (acmod_fe_mismatch(acmod, fe))
            goto error_out;
        fe_retain(fe);
        acmod->fe = fe;
    }
    else {
        /* Initialize a new front end. */
        acmod->fe = fe_init_auto_r(config);
        if (acmod->fe == NULL)
            goto error_out;
        if (acmod_fe_mismatch(acmod, acmod->fe))
            goto error_out;
    }
    if (fcb) {
        if (acmod_feat_mismatch(acmod, fcb))
            goto error_out;
        feat_retain(fcb);
        acmod->fcb = fcb;
    }
    else {
        /* Initialize a new fcb. */
        if (acmod_init_feat(acmod) < 0)
            goto error_out;
    }

    /* Load acoustic model parameters. */
    if (acmod_init_am(acmod) < 0)
        goto error_out;


    /* The MFCC buffer needs to be at least as large as the dynamic
     * feature window.  */
    acmod->n_mfc_alloc = acmod->fcb->window_size * 2 + 1;
    acmod->mfc_buf = (mfcc_t **)
        ckd_calloc_2d(acmod->n_mfc_alloc, acmod->fcb->cepsize,
                      sizeof(**acmod->mfc_buf));

    /* Feature buffer has to be at least as large as MFCC buffer. */
    acmod->n_feat_alloc = acmod->n_mfc_alloc + cmd_ln_int32_r(config, "-pl_window");
    acmod->feat_buf = feat_array_alloc(acmod->fcb, acmod->n_feat_alloc);
    acmod->framepos = ckd_calloc(acmod->n_feat_alloc, sizeof(*acmod->framepos));

    /* Senone computation stuff. */
    acmod->senone_scores = ckd_calloc(bin_mdef_n_sen(acmod->mdef),
                                                     sizeof(*acmod->senone_scores));
    acmod->senone_active_vec = bitvec_alloc(bin_mdef_n_sen(acmod->mdef));
    acmod->senone_active = ckd_calloc(bin_mdef_n_sen(acmod->mdef),
                                                     sizeof(*acmod->senone_active));
    acmod->log_zero = logmath_get_zero(acmod->lmath);
    acmod->compallsen = cmd_ln_boolean_r(config, "-compallsen");
    return acmod;

error_out:
    acmod_free(acmod);
    return NULL;
}

void
acmod_free(acmod_t *acmod)
{
    if (acmod == NULL)
        return;

    feat_free(acmod->fcb);
    fe_free(acmod->fe);
    cmd_ln_free_r(acmod->config);

    if (acmod->mfc_buf)
        ckd_free_2d((void **)acmod->mfc_buf);
    if (acmod->feat_buf)
        feat_array_free(acmod->feat_buf);

    if (acmod->mfcfh)
        fclose(acmod->mfcfh);
    if (acmod->rawfh)
        fclose(acmod->rawfh);
    if (acmod->senfh)
        fclose(acmod->senfh);

    ckd_free(acmod->framepos);
    ckd_free(acmod->senone_scores);
    ckd_free(acmod->senone_active_vec);
    ckd_free(acmod->senone_active);

    if (acmod->mdef)
        bin_mdef_free(acmod->mdef);
    if (acmod->tmat)
        tmat_free(acmod->tmat);
    if (acmod->mgau)
        ps_mgau_free(acmod->mgau);
    if (acmod->mllr)
        ps_mllr_free(acmod->mllr);
    
    ckd_free(acmod);
}

ps_mllr_t *
acmod_update_mllr(acmod_t *acmod, ps_mllr_t *mllr)
{
    if (acmod->mllr)
        ps_mllr_free(acmod->mllr);
    acmod->mllr = mllr;
    ps_mgau_transform(acmod->mgau, mllr);

    return mllr;
}

int
acmod_write_senfh_header(acmod_t *acmod, FILE *logfh)
{
    char nsenstr[64], logbasestr[64];

    sprintf(nsenstr, "%d", bin_mdef_n_sen(acmod->mdef));
    sprintf(logbasestr, "%f", logmath_get_base(acmod->lmath));
    return bio_writehdr(logfh,
                        "version", "0.1",
                        "mdef_file", cmd_ln_str_r(acmod->config, "-mdef"),
                        "n_sen", nsenstr,
                        "logbase", logbasestr, NULL);
}

int
acmod_set_senfh(acmod_t *acmod, FILE *logfh)
{
    if (acmod->senfh)
        fclose(acmod->senfh);
    acmod->senfh = logfh;
    if (logfh == NULL)
        return 0;
    return acmod_write_senfh_header(acmod, logfh);
}

int
acmod_set_mfcfh(acmod_t *acmod, FILE *logfh)
{
    int rv = 0;

    if (acmod->mfcfh)
        fclose(acmod->mfcfh);
    acmod->mfcfh = logfh;
    fwrite(&rv, 4, 1, acmod->mfcfh);
    return rv;
}

int
acmod_set_rawfh(acmod_t *acmod, FILE *logfh)
{
    if (acmod->rawfh)
        fclose(acmod->rawfh);
    acmod->rawfh = logfh;
    return 0;
}

void
acmod_grow_feat_buf(acmod_t *acmod, int nfr)
{
    if (nfr > MAX_N_FRAMES)
        E_FATAL("Decoder can not process more than %d frames at once, requested %d\n", 
                MAX_N_FRAMES, nfr);

    acmod->feat_buf = feat_array_realloc(acmod->fcb, acmod->feat_buf, acmod->n_feat_alloc, nfr);
    acmod->framepos = ckd_realloc(acmod->framepos,
                                  nfr * sizeof(*acmod->framepos));
    acmod->n_feat_alloc = nfr;
}

int
acmod_set_grow(acmod_t *acmod, int grow_feat)
{
    int tmp = acmod->grow_feat;
    acmod->grow_feat = grow_feat;

    /* Expand feat_buf to a reasonable size to start with. */
    if (grow_feat && acmod->n_feat_alloc < 128)
        acmod_grow_feat_buf(acmod, 128);

    return tmp;
}

int
acmod_start_utt(acmod_t *acmod)
{
    fe_start_utt(acmod->fe);
    acmod->state = ACMOD_STARTED;
    acmod->n_mfc_frame = 0;
    acmod->n_feat_frame = 0;
    acmod->mfc_outidx = 0;
    acmod->feat_outidx = 0;
    acmod->output_frame = 0;
    acmod->senscr_frame = -1;
    acmod->n_senone_active = 0;
    acmod->mgau->frame_idx = 0;
    return 0;
}

int
acmod_end_utt(acmod_t *acmod)
{
    int32 nfr = 0;

    acmod->state = ACMOD_ENDED;
    if (acmod->n_mfc_frame < acmod->n_mfc_alloc) {
        int inptr;
        /* Where to start writing them (circular buffer) */
        inptr = (acmod->mfc_outidx + acmod->n_mfc_frame) % acmod->n_mfc_alloc;
        /* nfr is always either zero or one. */
        fe_end_utt(acmod->fe, acmod->mfc_buf[inptr], &nfr);
        acmod->n_mfc_frame += nfr;
        /* Process whatever's left, and any leadout. */
        if (nfr)
            nfr = acmod_process_mfcbuf(acmod);
    }
    if (acmod->mfcfh) {
        int32 outlen, rv;
        outlen = (ftell(acmod->mfcfh) - 4) / 4;
        if (!WORDS_BIGENDIAN)
            SWAP_INT32(&outlen);
        /* Try to seek and write */
        if ((rv = fseek(acmod->mfcfh, 0, SEEK_SET)) == 0) {
            fwrite(&outlen, 4, 1, acmod->mfcfh);
        }
        fclose(acmod->mfcfh);
        acmod->mfcfh = NULL;
    }
    if (acmod->rawfh) {
        fclose(acmod->rawfh);
        acmod->rawfh = NULL;
    }

    if (acmod->senfh) {
        fclose(acmod->senfh);
        acmod->senfh = NULL;
    }

    return nfr;
}

static int
acmod_log_mfc(acmod_t *acmod,
              mfcc_t **cep, int n_frames)
{
    int i, n;
    int32 *ptr = (int32 *)cep[0];

    n = n_frames * feat_cepsize(acmod->fcb);
    /* Swap bytes. */
    if (!WORDS_BIGENDIAN) {
        for (i = 0; i < (n * sizeof(mfcc_t)); ++i) {
            SWAP_INT32(ptr + i);
        }
    }
    /* Write features. */
    if (fwrite(cep[0], sizeof(mfcc_t), n, acmod->mfcfh) != n) {
        E_ERROR_SYSTEM("Failed to write %d values to log file", n);
    }

    /* Swap them back. */
    if (!WORDS_BIGENDIAN) {
        for (i = 0; i < (n * sizeof(mfcc_t)); ++i) {
            SWAP_INT32(ptr + i);
        }
    }
    return 0;
}

static int
acmod_process_full_cep(acmod_t *acmod,
                       mfcc_t ***inout_cep,
                       int *inout_n_frames)
{
    int32 nfr;

    /* Write to log file. */
    if (acmod->mfcfh)
        acmod_log_mfc(acmod, *inout_cep, *inout_n_frames);

    /* Resize feat_buf to fit. */
    if (acmod->n_feat_alloc < *inout_n_frames) {
	    
	if (*inout_n_frames > MAX_N_FRAMES)
	    E_FATAL("Batch processing can not process more than %d frames at once, requested %d\n", 
		    MAX_N_FRAMES, *inout_n_frames);
    
        feat_array_free(acmod->feat_buf);
        acmod->feat_buf = feat_array_alloc(acmod->fcb, *inout_n_frames);
        acmod->n_feat_alloc = *inout_n_frames;
        acmod->n_feat_frame = 0;
        acmod->feat_outidx = 0;
    }
    /* Make dynamic features. */
    nfr = feat_s2mfc2feat_live(acmod->fcb, *inout_cep, inout_n_frames,
                               TRUE, TRUE, acmod->feat_buf);
    acmod->n_feat_frame = nfr;
    assert(acmod->n_feat_frame <= acmod->n_feat_alloc);
    *inout_cep += *inout_n_frames;
    *inout_n_frames = 0;
    return nfr;
}

static int
acmod_process_full_raw(acmod_t *acmod,
                       int16 const **inout_raw,
                       size_t *inout_n_samps)
{
    int32 nfr, ntail;
    mfcc_t **cepptr;

    /* Write to logging file if any. */
    if (acmod->rawfh)
        fwrite(*inout_raw, 2, *inout_n_samps, acmod->rawfh);
    /* Resize mfc_buf to fit. */
    if (fe_process_frames(acmod->fe, NULL, inout_n_samps, NULL, &nfr, NULL) < 0)
        return -1;
    if (acmod->n_mfc_alloc < nfr + 1) {
        ckd_free_2d(acmod->mfc_buf);
        acmod->mfc_buf = ckd_calloc_2d(nfr + 1, fe_get_output_size(acmod->fe),
                                       sizeof(**acmod->mfc_buf));
        acmod->n_mfc_alloc = nfr + 1;
    }
    acmod->n_mfc_frame = 0;
    acmod->mfc_outidx = 0;
    fe_start_utt(acmod->fe);
    if (fe_process_frames(acmod->fe, inout_raw, inout_n_samps, acmod->mfc_buf, &nfr, NULL) < 0)
        return -1;
    fe_end_utt(acmod->fe, acmod->mfc_buf[nfr], &ntail);
    nfr += ntail;

    cepptr = acmod->mfc_buf;
    nfr = acmod_process_full_cep(acmod, &cepptr, &nfr);
    acmod->n_mfc_frame = 0;
    return nfr;
}

/**
 * Process MFCCs that are in the internal buffer into features.
 */
static int32
acmod_process_mfcbuf(acmod_t *acmod)
{
    mfcc_t **mfcptr;
    int32 ncep;

    ncep = acmod->n_mfc_frame;
    /* Also do this in two parts because of the circular mfc_buf. */
    if (acmod->mfc_outidx + ncep > acmod->n_mfc_alloc) {
        int32 ncep1 = acmod->n_mfc_alloc - acmod->mfc_outidx;
        int saved_state = acmod->state;

        /* Make sure we don't end the utterance here. */
        if (acmod->state == ACMOD_ENDED)
            acmod->state = ACMOD_PROCESSING;
        mfcptr = acmod->mfc_buf + acmod->mfc_outidx;
        ncep1 = acmod_process_cep(acmod, &mfcptr, &ncep1, FALSE);
        /* It's possible that not all available frames were filled. */
        ncep -= ncep1;
        acmod->n_mfc_frame -= ncep1;
        acmod->mfc_outidx += ncep1;
        acmod->mfc_outidx %= acmod->n_mfc_alloc;
        /* Restore original state (could this really be the end) */
        acmod->state = saved_state;
    }
    mfcptr = acmod->mfc_buf + acmod->mfc_outidx;
    ncep = acmod_process_cep(acmod, &mfcptr, &ncep, FALSE);
    acmod->n_mfc_frame -= ncep;
    acmod->mfc_outidx += ncep;
    acmod->mfc_outidx %= acmod->n_mfc_alloc;
    return ncep;
}

int
acmod_process_raw(acmod_t *acmod,
		  int16 const **inout_raw,
		  size_t *inout_n_samps,
		  int full_utt)
{
    int32 ncep;

    /* If this is a full utterance, process it all at once. */
    if (full_utt)
        return acmod_process_full_raw(acmod, inout_raw, inout_n_samps);

    /* Append MFCCs to the end of any that are previously in there
     * (in practice, there will probably be none) */
    if (inout_n_samps && *inout_n_samps) {
        int16 const *prev_audio_inptr = *inout_raw;
        int inptr;

        /* Total number of frames available. */
        ncep = acmod->n_mfc_alloc - acmod->n_mfc_frame;
        /* Where to start writing them (circular buffer) */
        inptr = (acmod->mfc_outidx + acmod->n_mfc_frame) % acmod->n_mfc_alloc;

        /* Write them in two (or more) parts if there is wraparound. */
        while (inptr + ncep > acmod->n_mfc_alloc) {
            int32 ncep1 = acmod->n_mfc_alloc - inptr;
            if (fe_process_frames(acmod->fe, inout_raw, inout_n_samps, acmod->mfc_buf + inptr, &ncep1, NULL) < 0)
                return -1;
            /* Write to logging file if any. */
            if (acmod->rawfh) {
                fwrite(prev_audio_inptr, 2,
                       *inout_raw - prev_audio_inptr,
                       acmod->rawfh);
                prev_audio_inptr = *inout_raw;
            }
            /* ncep1 now contains the number of frames actually
             * processed.  This is a good thing, but it means we
             * actually still might have some room left at the end of
             * the buffer, hence the while loop.  Unfortunately it
             * also means that in the case where we are really
             * actually done, we need to get out totally, hence the
             * goto. */
            acmod->n_mfc_frame += ncep1;
            ncep -= ncep1;
            inptr += ncep1;
            inptr %= acmod->n_mfc_alloc;
            if (ncep1 == 0)
                goto alldone;
        }
        assert(inptr + ncep <= acmod->n_mfc_alloc);
        if (fe_process_frames(acmod->fe, inout_raw, inout_n_samps, acmod->mfc_buf + inptr, &ncep, NULL) < 0)
            return -1;
        /* Write to logging file if any. */
        if (acmod->rawfh) {
            fwrite(prev_audio_inptr, 2,
                   *inout_raw - prev_audio_inptr, acmod->rawfh);
            prev_audio_inptr = *inout_raw;
        }
        acmod->n_mfc_frame += ncep;
    alldone:
        ;
    }

    /* Hand things off to acmod_process_cep. */
    return acmod_process_mfcbuf(acmod);
}

int
acmod_process_cep(acmod_t *acmod,
		  mfcc_t ***inout_cep,
		  int *inout_n_frames,
		  int full_utt)
{
    int32 nfeat, ncep, inptr;
    int orig_n_frames;

    /* If this is a full utterance, process it all at once. */
    if (full_utt)
        return acmod_process_full_cep(acmod, inout_cep, inout_n_frames);

    /* Write to log file. */
    if (acmod->mfcfh)
        acmod_log_mfc(acmod, *inout_cep, *inout_n_frames);

    /* Maximum number of frames we're going to generate. */
    orig_n_frames = ncep = nfeat = *inout_n_frames;

    /* FIXME: This behaviour isn't guaranteed... */
    if (acmod->state == ACMOD_ENDED)
        nfeat += feat_window_size(acmod->fcb);
    else if (acmod->state == ACMOD_STARTED)
        nfeat -= feat_window_size(acmod->fcb);

    /* Clamp number of features to fit available space. */
    if (nfeat > acmod->n_feat_alloc - acmod->n_feat_frame) {
        /* Grow it as needed - we have to grow it at the end of an
         * utterance because we can't return a short read there. */
        if (acmod->grow_feat || acmod->state == ACMOD_ENDED)
            acmod_grow_feat_buf(acmod, acmod->n_feat_alloc + nfeat);
        else
            ncep -= (nfeat - (acmod->n_feat_alloc - acmod->n_feat_frame));
    }

    /* Where to start writing in the feature buffer. */
    if (acmod->grow_feat) {
        /* Grow to avoid wraparound if grow_feat == TRUE. */
        inptr = acmod->feat_outidx + acmod->n_feat_frame;
        while (inptr + nfeat >= acmod->n_feat_alloc)
            acmod_grow_feat_buf(acmod, acmod->n_feat_alloc * 2);
    }
    else {
        inptr = (acmod->feat_outidx + acmod->n_feat_frame) % acmod->n_feat_alloc;
    }


    /* FIXME: we can't split the last frame drop properly to be on the bounary, so just return */
    if (inptr + nfeat > acmod->n_feat_alloc && acmod->state == ACMOD_ENDED) {
	*inout_n_frames -= ncep;
	*inout_cep += ncep;
	return 0;
    }

    /* Write them in two parts if there is wraparound. */
    if (inptr + nfeat > acmod->n_feat_alloc) {
        int32 ncep1 = acmod->n_feat_alloc - inptr;

        /* Make sure we don't end the utterance here. */
        nfeat = feat_s2mfc2feat_live(acmod->fcb, *inout_cep,
                                     &ncep1,
                                     (acmod->state == ACMOD_STARTED),
                                     FALSE,
                                     acmod->feat_buf + inptr);
        if (nfeat < 0)
            return -1;
        /* Move the output feature pointer forward. */
        acmod->n_feat_frame += nfeat;
        assert(acmod->n_feat_frame <= acmod->n_feat_alloc);
        inptr += nfeat;
        inptr %= acmod->n_feat_alloc;
        /* Move the input feature pointers forward. */
        *inout_n_frames -= ncep1;
        *inout_cep += ncep1;
        ncep -= ncep1;
    }

    nfeat = feat_s2mfc2feat_live(acmod->fcb, *inout_cep,
                                 &ncep,
                                 (acmod->state == ACMOD_STARTED),
                                 (acmod->state == ACMOD_ENDED),
                                 acmod->feat_buf + inptr);
    if (nfeat < 0)
        return -1;
    acmod->n_feat_frame += nfeat;
    assert(acmod->n_feat_frame <= acmod->n_feat_alloc);
    /* Move the input feature pointers forward. */
    *inout_n_frames -= ncep;
    *inout_cep += ncep;
    if (acmod->state == ACMOD_STARTED)
        acmod->state = ACMOD_PROCESSING;
    return orig_n_frames - *inout_n_frames;
}

int
acmod_process_feat(acmod_t *acmod,
		   mfcc_t **feat)
{
    int i, inptr;

    if (acmod->n_feat_frame == acmod->n_feat_alloc) {
        if (acmod->grow_feat)
            acmod_grow_feat_buf(acmod, acmod->n_feat_alloc * 2);
        else
            return 0;
    }

    if (acmod->grow_feat) {
        /* Grow to avoid wraparound if grow_feat == TRUE. */
        inptr = acmod->feat_outidx + acmod->n_feat_frame;
        while (inptr + 1 >= acmod->n_feat_alloc)
            acmod_grow_feat_buf(acmod, acmod->n_feat_alloc * 2);
    }
    else {
        inptr = (acmod->feat_outidx + acmod->n_feat_frame) % acmod->n_feat_alloc;
    }
    for (i = 0; i < feat_dimension1(acmod->fcb); ++i)
        memcpy(acmod->feat_buf[inptr][i],
               feat[i], feat_dimension2(acmod->fcb, i) * sizeof(**feat));
    ++acmod->n_feat_frame;
    assert(acmod->n_feat_frame <= acmod->n_feat_alloc);

    return 1;
}

static int
acmod_read_senfh_header(acmod_t *acmod)
{
    char **name, **val;
    int32 swap;
    int i;

    if (bio_readhdr(acmod->insenfh, &name, &val, &swap) < 0)
        goto error_out;
    for (i = 0; name[i] != NULL; ++i) {
        if (!strcmp(name[i], "n_sen")) {
            if (atoi(val[i]) != bin_mdef_n_sen(acmod->mdef)) {
                E_ERROR("Number of senones in senone file (%d) does not match mdef (%d)\n",
                        atoi(val[i]), bin_mdef_n_sen(acmod->mdef));
                goto error_out;
            }
        }
        if (!strcmp(name[i], "logbase")) {
            if (abs(atof(val[i]) - logmath_get_base(acmod->lmath)) > 0.001) {
                E_ERROR("Logbase in senone file (%f) does not match acmod (%f)\n",
                        atof(val[i]), logmath_get_base(acmod->lmath));
                goto error_out;
            }
        }
    }
    acmod->insen_swap = swap;
    bio_hdrarg_free(name, val);
    return 0;
error_out:
    bio_hdrarg_free(name, val);
    return -1;
}

int
acmod_set_insenfh(acmod_t *acmod, FILE *senfh)
{
    acmod->insenfh = senfh;
    if (senfh == NULL) {
        acmod->n_feat_frame = 0;
        acmod->compallsen = cmd_ln_boolean_r(acmod->config, "-compallsen");
        return 0;
    }
    acmod->compallsen = TRUE;
    return acmod_read_senfh_header(acmod);
}

int
acmod_rewind(acmod_t *acmod)
{
    /* If the feature buffer is circular, this is not possible. */
    if (acmod->output_frame > acmod->n_feat_alloc) {
        E_ERROR("Circular feature buffer cannot be rewound (output frame %d, alloc %d)\n",
               acmod->output_frame, acmod->n_feat_alloc);
        return -1;
    }

    /* Frames consumed + frames available */
    acmod->n_feat_frame = acmod->output_frame + acmod->n_feat_frame;

    /* Reset output pointers. */
    acmod->feat_outidx = 0;
    acmod->output_frame = 0;
    acmod->senscr_frame = -1;
    acmod->mgau->frame_idx = 0;

    return 0;
}

int
acmod_advance(acmod_t *acmod)
{
    /* Advance the output pointers. */
    if (++acmod->feat_outidx == acmod->n_feat_alloc)
        acmod->feat_outidx = 0;
    --acmod->n_feat_frame;
    ++acmod->mgau->frame_idx;

    return ++acmod->output_frame;
}

int
acmod_write_scores(acmod_t *acmod, int n_active, uint8 const *active,
                   int16 const *senscr, FILE *senfh)
{
    int16 n_active2;

    /* Uncompressed frame format:
     *
     * (2 bytes) n_active: Number of active senones
     * If all senones active:
     * (n_active * 2 bytes) scores of active senones
     *
     * Otherwise:
     * (2 bytes) n_active: Number of active senones
     * (n_active bytes) deltas to active senones
     * (n_active * 2 bytes) scores of active senones
     */
    n_active2 = n_active;
    if (fwrite(&n_active2, 2, 1, senfh) != 1)
        goto error_out;
    if (n_active == bin_mdef_n_sen(acmod->mdef)) {
        if (fwrite(senscr, 2, n_active, senfh) != n_active)
            goto error_out;
    }
    else {
        int i, n;
        if (fwrite(active, 1, n_active, senfh) != n_active)
            goto error_out;
        for (i = n = 0; i < n_active; ++i) {
            n += active[i];
            if (fwrite(senscr + n, 2, 1, senfh) != 1)
                goto error_out;
        }
    }
    return 0;
error_out:
    E_ERROR_SYSTEM("Failed to write frame to senone file");
    return -1;
}

/**
 * Internal version, used for reading previous frames in acmod_score()
 */
static int
acmod_read_scores_internal(acmod_t *acmod)
{
    FILE *senfh = acmod->insenfh;
    int16 n_active;
    int rv;

    if (acmod->n_feat_frame == acmod->n_feat_alloc) {
        if (acmod->grow_feat)
            acmod_grow_feat_buf(acmod, acmod->n_feat_alloc * 2);
        else
            return 0;
    }

    if (senfh == NULL)
        return -1;
    if ((rv = fread(&n_active, 2, 1, senfh)) < 0)
        goto error_out;
    else if (rv == 0)
        return 0;

    acmod->n_senone_active = n_active;
    if (acmod->n_senone_active == bin_mdef_n_sen(acmod->mdef)) {
        if ((rv = fread(acmod->senone_scores, 2,
                        acmod->n_senone_active, senfh)) < 0)
            goto error_out;
        else if (rv != acmod->n_senone_active)
            return 0;
    }
    else {
        int i, n;
        if ((rv = fread(acmod->senone_active, 1,
                        acmod->n_senone_active, senfh)) < 0)
            goto error_out;
        else if (rv != acmod->n_senone_active)
            return 0;
        for (i = 0, n = 0; i < acmod->n_senone_active; ++i) {
            int j, sen = n + acmod->senone_active[i];
            for (j = n + 1; j < sen; ++j)
                acmod->senone_scores[j] = SENSCR_DUMMY;
            if ((rv = fread(acmod->senone_scores + sen, 2, 1, senfh)) < 0)
                goto error_out;
            else if (rv == 0)
                return 0;
            n = sen;
        }
        ++n;
        while (n < bin_mdef_n_sen(acmod->mdef))
            acmod->senone_scores[n++] = SENSCR_DUMMY;
    }
    return 1;
error_out:
    E_ERROR_SYSTEM("Failed to read frame from senone file");
    return -1;
}

int
acmod_read_scores(acmod_t *acmod)
{
    int inptr, rv;

    if (acmod->grow_feat) {
        /* Grow to avoid wraparound if grow_feat == TRUE. */
        inptr = acmod->feat_outidx + acmod->n_feat_frame;
        /* Has to be +1, otherwise, next time acmod_advance() is
         * called, this will wrap around. */
        while (inptr + 1 >= acmod->n_feat_alloc)
            acmod_grow_feat_buf(acmod, acmod->n_feat_alloc * 2);
    }
    else {
        inptr = (acmod->feat_outidx + acmod->n_feat_frame) % acmod->n_feat_alloc;
    }

    if ((rv = acmod_read_scores_internal(acmod)) != 1)
        return rv;

    /* Set acmod->senscr_frame appropriately so that these scores
       get reused below in acmod_score(). */
    acmod->senscr_frame = acmod->output_frame + acmod->n_feat_frame;

    E_DEBUG(1,("Frame %d has %d active states\n",
               acmod->senscr_frame, acmod->n_senone_active));

    /* Increment the "feature frame counter" and record the file
     * position for the relevant frame in the (possibly circular)
     * buffer. */
    ++acmod->n_feat_frame;
    acmod->framepos[inptr] = ftell(acmod->insenfh);

    return 1;
}

static int
calc_frame_idx(acmod_t *acmod, int *inout_frame_idx)
{
    int frame_idx;

    /* Calculate the absolute frame index to be scored. */
    if (inout_frame_idx == NULL)
        frame_idx = acmod->output_frame;
    else if (*inout_frame_idx < 0)
        frame_idx = acmod->output_frame + 1 + *inout_frame_idx;
    else
        frame_idx = *inout_frame_idx;

    return frame_idx;
}

static int
calc_feat_idx(acmod_t *acmod, int frame_idx)
{
    int n_backfr, feat_idx;

    n_backfr = acmod->n_feat_alloc - acmod->n_feat_frame;
    if (frame_idx < 0 || acmod->output_frame - frame_idx > n_backfr) {
        E_ERROR("Frame %d outside queue of %d frames, %d alloc (%d > %d), cannot score\n",
                frame_idx, acmod->n_feat_frame, acmod->n_feat_alloc,
                acmod->output_frame - frame_idx, n_backfr);
        return -1;
    }

    /* Get the index in feat_buf/framepos of the frame to be scored. */
    feat_idx = ((acmod->feat_outidx + frame_idx - acmod->output_frame)
                % acmod->n_feat_alloc);
    if (feat_idx < 0) feat_idx += acmod->n_feat_alloc;

    return feat_idx;
}

mfcc_t **
acmod_get_frame(acmod_t *acmod, int *inout_frame_idx)
{
    int frame_idx, feat_idx;

    /* Calculate the absolute frame index requested. */
    frame_idx = calc_frame_idx(acmod, inout_frame_idx);

    /* Calculate position of requested frame in circular buffer. */
    if ((feat_idx = calc_feat_idx(acmod, frame_idx)) < 0)
        return NULL;

    if (inout_frame_idx)
        *inout_frame_idx = frame_idx;

    return acmod->feat_buf[feat_idx];
}

int16 const *
acmod_score(acmod_t *acmod, int *inout_frame_idx)
{
    int frame_idx, feat_idx;

    /* Calculate the absolute frame index to be scored. */
    frame_idx = calc_frame_idx(acmod, inout_frame_idx);

    /* If all senones are being computed, or we are using a senone file,
       then we can reuse existing scores. */
    if ((acmod->compallsen || acmod->insenfh)
        && frame_idx == acmod->senscr_frame) {
        if (inout_frame_idx)
            *inout_frame_idx = frame_idx;
        return acmod->senone_scores;
    }

    /* Calculate position of requested frame in circular buffer. */
    if ((feat_idx = calc_feat_idx(acmod, frame_idx)) < 0)
        return NULL;

    /* If there is an input senone file locate the appropriate frame and read it. */
    if (acmod->insenfh) {
        fseek(acmod->insenfh, acmod->framepos[feat_idx], SEEK_SET);
        if (acmod_read_scores_internal(acmod) < 0)
            return NULL;
    }
    else {
        /* Build active senone list. */
        acmod_flags2list(acmod);

        /* Generate scores for the next available frame */
        ps_mgau_frame_eval(acmod->mgau,
                           acmod->senone_scores,
                           acmod->senone_active,
                           acmod->n_senone_active,
                           acmod->feat_buf[feat_idx],
                           frame_idx,
                           acmod->compallsen);
    }

    if (inout_frame_idx)
        *inout_frame_idx = frame_idx;
    acmod->senscr_frame = frame_idx;

    /* Dump scores to the senone dump file if one exists. */
    if (acmod->senfh) {
        if (acmod_write_scores(acmod, acmod->n_senone_active,
                               acmod->senone_active,
                               acmod->senone_scores,
                               acmod->senfh) < 0)
            return NULL;
        E_DEBUG(1,("Frame %d has %d active states\n", frame_idx, acmod->n_senone_active));
    }

    return acmod->senone_scores;
}

int
acmod_best_score(acmod_t *acmod, int *out_best_senid)
{
    int i, best;

    best = SENSCR_DUMMY;
    if (acmod->compallsen) {
        for (i = 0; i < bin_mdef_n_sen(acmod->mdef); ++i) {
            if (acmod->senone_scores[i] < best) {
                best = acmod->senone_scores[i];
                *out_best_senid = i;
            }
        }
    }
    else {
        int16 *senscr;
        senscr = acmod->senone_scores;
        for (i = 0; i < acmod->n_senone_active; ++i) {
            senscr += acmod->senone_active[i];
            if (*senscr < best) {
                best = *senscr;
                *out_best_senid = i;
            }
        }
    }
    return best;
}


void
acmod_clear_active(acmod_t *acmod)
{
    if (acmod->compallsen)
        return;
    bitvec_clear_all(acmod->senone_active_vec, bin_mdef_n_sen(acmod->mdef));
    acmod->n_senone_active = 0;
}

#define MPX_BITVEC_SET(a,h,i)                                   \
    if (hmm_mpx_ssid(h,i) != BAD_SSID)                          \
        bitvec_set((a)->senone_active_vec, hmm_mpx_senid(h,i))
#define NONMPX_BITVEC_SET(a,h,i)                                        \
    bitvec_set((a)->senone_active_vec,                                  \
               hmm_nonmpx_senid(h,i))

void
acmod_activate_hmm(acmod_t *acmod, hmm_t *hmm)
{
    int i;

    if (acmod->compallsen)
        return;
    if (hmm_is_mpx(hmm)) {
        switch (hmm_n_emit_state(hmm)) {
        case 5:
            MPX_BITVEC_SET(acmod, hmm, 4);
            MPX_BITVEC_SET(acmod, hmm, 3);
        case 3:
            MPX_BITVEC_SET(acmod, hmm, 2);
            MPX_BITVEC_SET(acmod, hmm, 1);
            MPX_BITVEC_SET(acmod, hmm, 0);
            break;
        default:
            for (i = 0; i < hmm_n_emit_state(hmm); ++i) {
                MPX_BITVEC_SET(acmod, hmm, i);
            }
        }
    }
    else {
        switch (hmm_n_emit_state(hmm)) {
        case 5:
            NONMPX_BITVEC_SET(acmod, hmm, 4);
            NONMPX_BITVEC_SET(acmod, hmm, 3);
        case 3:
            NONMPX_BITVEC_SET(acmod, hmm, 2);
            NONMPX_BITVEC_SET(acmod, hmm, 1);
            NONMPX_BITVEC_SET(acmod, hmm, 0);
            break;
        default:
            for (i = 0; i < hmm_n_emit_state(hmm); ++i) {
                NONMPX_BITVEC_SET(acmod, hmm, i);
            }
        }
    }
}

int32
acmod_flags2list(acmod_t *acmod)
{
    int32 w, l, n, b, total_dists, total_words, extra_bits;
    bitvec_t *flagptr;

    total_dists = bin_mdef_n_sen(acmod->mdef);
    if (acmod->compallsen) {
        acmod->n_senone_active = total_dists;
        return total_dists;
    }
    total_words = total_dists / BITVEC_BITS;
    extra_bits = total_dists % BITVEC_BITS;
    w = n = l = 0;
    for (flagptr = acmod->senone_active_vec; w < total_words; ++w, ++flagptr) {
        if (*flagptr == 0)
            continue;
        for (b = 0; b < BITVEC_BITS; ++b) {
            if (*flagptr & (1UL << b)) {
                int32 sen = w * BITVEC_BITS + b;
                int32 delta = sen - l;
                /* Handle excessive deltas "lossily" by adding a few
                   extra senones to bridge the gap. */
                while (delta > 255) {
                    acmod->senone_active[n++] = 255;
                    delta -= 255;
                }
                acmod->senone_active[n++] = delta;
                l = sen;
            }
        }
    }

    for (b = 0; b < extra_bits; ++b) {
        if (*flagptr & (1UL << b)) {
            int32 sen = w * BITVEC_BITS + b;
            int32 delta = sen - l;
            /* Handle excessive deltas "lossily" by adding a few
               extra senones to bridge the gap. */
            while (delta > 255) {
                acmod->senone_active[n++] = 255;
                delta -= 255;
            }
            acmod->senone_active[n++] = delta;
            l = sen;
        }
    }

    acmod->n_senone_active = n;
    E_DEBUG(1, ("acmod_flags2list: %d active in frame %d\n",
                acmod->n_senone_active, acmod->output_frame));
    return n;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file acmod.h Acoustic model structures for PocketSphinx.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __ACMOD_H__
#define __ACMOD_H__

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/bitvec.h>
#include <sphinxbase/err.h>
#include <sphinxbase/prim_type.h>

/* Local headers. */
#include "ps_mllr.h"
#include "bin_mdef.h"
#include "tmat.h"
#include "hmm.h"

/**
 * States in utterance processing.
 */
typedef enum acmod_state_e {
    ACMOD_IDLE,		/**< Not in an utterance. */
    ACMOD_STARTED,      /**< Utterance started, no data yet. */
    ACMOD_PROCESSING,   /**< Utterance in progress. */
    ACMOD_ENDED         /**< Utterance ended, still buffering. */
} acmod_state_t;

/**
 * Dummy senone score value for unintentionally active states.
 */
#define SENSCR_DUMMY 0x7fff

/**
 * Feature space linear transform structure.
 */
struct ps_mllr_s {
    int refcnt;     /**< Reference count. */
    int n_class;    /**< Number of MLLR classes. */
    int n_feat;     /**< Number of feature streams. */
    int *veclen;    /**< Length of input vectors for each stream. */
    float32 ****A;  /**< Rotation part of mean transformations. */
    float32 ***b;   /**< Bias part of mean transformations. */
    float32 ***h;   /**< Diagonal transformation of variances. */
    int32 *cb2mllr; /**< Mapping from codebooks to transformations. */
};

/**
 * Acoustic model parameter structure. 
 */
typedef struct ps_mgau_s ps_mgau_t;

typedef struct ps_mgaufuncs_s {
    char const *name;

    int (*frame_eval)(ps_mgau_t *mgau,
                      int16 *senscr,
                      uint8 *senone_active,
                      int32 n_senone_active,
                      mfcc_t ** feat,
                      int32 frame,
                      int32 compallsen);
    int (*transform)(ps_mgau_t *mgau,
                     ps_mllr_t *mllr);
    void (*free)(ps_mgau_t *mgau);
} ps_mgaufuncs_t;    

struct ps_mgau_s {
    ps_mgaufuncs_t *vt;  /**< vtable of mgau functions. */
    int frame_idx;       /**< frame counter. */
};

#define ps_mgau_base(mg) ((ps_mgau_t *)(mg))
#define ps_mgau_frame_eval(mg,senscr,senone_active,n_senone_active,feat,frame,compallsen) \
    (*ps_mgau_base(mg)->vt->frame_eval)                                 \
    (mg, senscr, senone_active, n_senone_active, feat, frame, compallsen)
#define ps_mgau_transform(mg, mllr)                                  \
    (*ps_mgau_base(mg)->vt->transform)(mg, mllr)
#define ps_mgau_free(mg)                                  \
    (*ps_mgau_base(mg)->vt->free)(mg)

/**
 * Acoustic model structure.
 *
 * This object encapsulates all stages of acoustic processing, from
 * raw audio input to acoustic score output.  The reason for grouping
 * all of these modules together is that they all have to "agree" in
 * their parameterizations, and the configuration of the acoustic and
 * dynamic feature computation is completely dependent on the
 * parameters used to build the original acoustic model (which should
 * by now always be specified in a feat.params file).
 *
 * Because there is not a one-to-one correspondence from blocks of
 * input audio or frames of input features to frames of acoustic
 * scores (due to dynamic feature calculation), results may not be
 * immediately available after input, and the output results will not
 * correspond to the last piece of data input.
 *
 * TODO: In addition, this structure serves the purpose of queueing
 * frames of features (and potentially also scores in the future) for
 * asynchronous passes of recognition operating in parallel.
 */
struct acmod_s {
    /* Global objects, not retained. */
    cmd_ln_t *config;          /**< Configuration. */
    logmath_t *lmath;          /**< Log-math computation. */
    glist_t strings;           /**< Temporary acoustic model filenames. */

    /* Feature computation: */
    fe_t *fe;                  /**< Acoustic feature computation. */
    feat_t *fcb;               /**< Dynamic feature computation. */

    /* Model parameters: */
    bin_mdef_t *mdef;          /**< Model definition. */
    tmat_t *tmat;              /**< Transition matrices. */
    ps_mgau_t *mgau;           /**< Model parameters. */
    ps_mllr_t *mllr;           /**< Speaker transformation. */

    /* Senone scoring: */
    int16 *senone_scores;      /**< GMM scores for current frame. */
    bitvec_t *senone_active_vec; /**< Active GMMs in current frame. */
    uint8 *senone_active;      /**< Array of deltas to active GMMs. */
    int senscr_frame;          /**< Frame index for senone_scores. */
    int n_senone_active;       /**< Number of active GMMs. */
    int log_zero;              /**< Zero log-probability value. */

    /* Utterance processing: */
    mfcc_t **mfc_buf;   /**< Temporary buffer of acoustic features. */
    mfcc_t ***feat_buf; /**< Temporary buffer of dynamic features. */
    FILE *rawfh;        /**< File for writing raw audio data. */
    FILE *mfcfh;        /**< File for writing acoustic feature data. */
    FILE *senfh;        /**< File for writing senone score data. */
    FILE *insenfh;	/**< Input senone score file. */
    long *framepos;     /**< File positions of recent frames in senone file. */

    /* A whole bunch of flags and counters: */
    uint8 state;        /**< State of utterance processing. */
    uint8 compallsen;   /**< Compute all senones? */
    uint8 grow_feat;    /**< Whether to grow feat_buf. */
    uint8 insen_swap;   /**< Whether to swap input senone score. */

    frame_idx_t output_frame; /**< Index of next frame of dynamic features. */
    frame_idx_t n_mfc_alloc;  /**< Number of frames allocated in mfc_buf */
    frame_idx_t n_mfc_frame;  /**< Number of frames active in mfc_buf */
    frame_idx_t mfc_outidx;   /**< Start of active frames in mfc_buf */
    frame_idx_t n_feat_alloc; /**< Number of frames allocated in feat_buf */
    frame_idx_t n_feat_frame; /**< Number of frames active in feat_buf */
    frame_idx_t feat_outidx;  /**< Start of active frames in feat_buf */
};
typedef struct acmod_s acmod_t;

/**
 * Initialize an acoustic model.
 *
 * @param config a command-line object containing parameters.  This
 *               pointer is not retained by this object.
 * @param lmath global log-math parameters.
 * @param fe a previously-initialized acoustic feature module to use,
 *           or NULL to create one automatically.  If this is supplied
 *           and its parameters do not match those in the acoustic
 *           model, this function will fail.  This pointer is not retained.
 * @param fe a previously-initialized dynamic feature module to use,
 *           or NULL to create one automatically.  If this is supplied
 *           and its parameters do not match those in the acoustic
 *           model, this function will fail.  This pointer is not retained.
 * @return a newly initialized acmod_t, or NULL on failure.
 */
acmod_t *acmod_init(cmd_ln_t *config, logmath_t *lmath, fe_t *fe, feat_t *fcb);

/**
 * Adapt acoustic model using a linear transform.
 *
 * @param mllr The new transform to use, or NULL to update the existing
 *              transform.  The decoder retains ownership of this pointer,
 *              so you should not attempt to free it manually.  Use
 *              ps_mllr_retain() if you wish to reuse it
 *              elsewhere.
 * @return The updated transform object for this decoder, or
 *         NULL on failure.
 */
ps_mllr_t *acmod_update_mllr(acmod_t *acmod, ps_mllr_t *mllr);

/**
 * Start logging senone scores to a filehandle.
 *
 * @param acmod Acoustic model object.
 * @param logfh Filehandle to log to.
 * @return 0 for success, <0 on error.
 */
int acmod_set_senfh(acmod_t *acmod, FILE *senfh);

/**
 * Start logging MFCCs to a filehandle.
 *
 * @param acmod Acoustic model object.
 * @param logfh Filehandle to log to.
 * @return 0 for success, <0 on error.
 */
int acmod_set_mfcfh(acmod_t *acmod, FILE *logfh);

/**
 * Start logging raw audio to a filehandle.
 *
 * @param acmod Acoustic model object.
 * @param logfh Filehandle to log to.
 * @return 0 for success, <0 on error.
 */
int acmod_set_rawfh(acmod_t *acmod, FILE *logfh);

/**
 * Finalize an acoustic model.
 */
void acmod_free(acmod_t *acmod);

/**
 * Mark the start of an utterance.
 */
int acmod_start_utt(acmod_t *acmod);

/**
 * Mark the end of an utterance.
 */
int acmod_end_utt(acmod_t *acmod);

/**
 * Rewind the current utterance, allowing it to be rescored.
 *
 * After calling this function, the internal frame index is reset, and
 * acmod_score() will return scores starting at the first frame of the
 * current utterance.  Currently, acmod_set_grow() must have been
 * called to enable growing the feature buffer in order for this to
 * work.  In the future, senone scores may be cached instead.
 *
 * @return 0 for success, <0 for failure (if the utterance can't be
 *         rewound due to no feature or score data available)
 */
int acmod_rewind(acmod_t *acmod);

/**
 * Advance the frame index.
 *
 * This function moves to the next frame of input data.  Subsequent
 * calls to acmod_score() will return scores for that frame, until the
 * next call to acmod_advance().
 *
 * @return New frame index.
 */
int acmod_advance(acmod_t *acmod);

/**
 * Set memory allocation policy for utterance processing.
 *
 * @param grow_feat If non-zero, the internal dynamic feature buffer
 * will expand as necessary to encompass any amount of data fed to the
 * model.
 * @return previous allocation policy.
 */
int acmod_set_grow(acmod_t *acmod, int grow_feat);

/**
 * TODO: Set queue length for utterance processing.
 *
 * This function allows multiple concurrent passes of search to
 * operate on different parts of the utterance.
 */

/**
 * Feed raw audio data to the acoustic model for scoring.
 *
 * @param inout_raw In: Pointer to buffer of raw samples
 *                  Out: Pointer to next sample to be read
 * @param inout_n_samps In: Number of samples available
 *                      Out: Number of samples remaining
 * @param full_utt If non-zero, this block represents a full
 *                 utterance and should be processed as such.
 * @return Number of frames of data processed.
 */
int acmod_process_raw(acmod_t *acmod,
                      int16 const **inout_raw,
                      size_t *inout_n_samps,
                      int full_utt);

/**
 * Feed acoustic feature data into the acoustic model for scoring.
 *
 * @param inout_cep In: Pointer to buffer of features
 *                  Out: Pointer to next frame to be read
 * @param inout_n_frames In: Number of frames available
 *                      Out: Number of frames remaining
 * @param full_utt If non-zero, this block represents a full
 *                 utterance and should be processed as such.
 * @return Number of frames of data processed.
 */
int acmod_process_cep(acmod_t *acmod,
                      mfcc_t ***inout_cep,
                      int *inout_n_frames,
                      int full_utt);

/**
 * Feed dynamic feature data into the acoustic model for scoring.
 *
 * Unlike acmod_process_raw() and acmod_process_cep(), this function
 * accepts a single frame at a time.  This is because there is no need
 * to do buffering when using dynamic features as input.  However, if
 * the dynamic feature buffer is full, this function will fail, so you
 * should either always check the return value, or always pair a call
 * to it with a call to acmod_score().
 *
 * @param feat Pointer to one frame of dynamic features.
 * @return Number of frames processed (either 0 or 1).
 */
int acmod_process_feat(acmod_t *acmod,
                       mfcc_t **feat);

/**
 * Set up a senone score dump file for input.
 *
 * @param insenfh File handle of dump file
 * @return 0 for success, <0 for failure
 */
int acmod_set_insenfh(acmod_t *acmod, FILE *insenfh);

/**
 * Read one frame of scores from senone score dump file.
 *
 * @return Number of frames read or <0 on error.
 */
int acmod_read_scores(acmod_t *acmod);

/**
 * Get a frame of dynamic feature data.
 *
 * @param inout_frame_idx Input: frame index to get, or NULL
 *                        to obtain features for the most recent frame.
 *                        Output: frame index corresponding to this
 *                        set of features.
 * @return Feature array, or NULL if requested frame is not available.
 */
mfcc_t **acmod_get_frame(acmod_t *acmod, int *inout_frame_idx);

/**
 * Score one frame of data.
 *
 * @param inout_frame_idx Input: frame index to score, or NULL
 *                        to obtain scores for the most recent frame.
 *                        Output: frame index corresponding to this
 *                        set of scores.
 * @return Array of senone scores for this frame, or NULL if no frame
 *         is available for scoring (such as if a frame index is
 *         requested that is not yet or no longer available).  The
 *         data pointed to persists only until the next call to
 *         acmod_score() or acmod_advance().
 */
int16 const *acmod_score(acmod_t *acmod,
                         int *inout_frame_idx);

/**
 * Write senone dump file header.
 */
int acmod_write_senfh_header(acmod_t *acmod, FILE *logfh);

/**
 * Write a frame of senone scores to a dump file.
 */
int acmod_write_scores(acmod_t *acmod, int n_active, uint8 const *active,
                       int16 const *senscr, FILE *senfh);


/**
 * Get best score and senone index for current frame.
 */
int acmod_best_score(acmod_t *acmod, int *out_best_senid);

/**
 * Clear set of active senones.
 */
void acmod_clear_active(acmod_t *acmod);

/**
 * Activate senones associated with an HMM.
 */
void acmod_activate_hmm(acmod_t *acmod, hmm_t *hmm);

/**
 * Activate a single senone.
 */
#define acmod_activate_sen(acmod, sen) bitvec_set((acmod)->senone_active_vec, sen)

/**
 * Build active list from 
 */
int32 acmod_flags2list(acmod_t *acmod);

#endif /* __ACMOD_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2005 Carnegie Mellon University.  All rights 
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*********************************************************************
 *
 * File: bin_mdef.c
 * 
 * Description: 
 *	Binary format model definition files, with support for
 *	heterogeneous topologies and variable-size N-phones
 *
 * Author: 
 * 	David Huggins-Daines <dhuggins@cs.cmu.edu>
 *********************************************************************/

/* System headers. */
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/byteorder.h>
#include <sphinxbase/case.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "mdef.h"
#include "bin_mdef.h"

bin_mdef_t *
bin_mdef_read_text(cmd_ln_t *config, const char *filename)
{
    bin_mdef_t *bmdef;
    mdef_t *mdef;
    int i, nodes, ci_idx, lc_idx, rc_idx;
    int nchars;

    if ((mdef = mdef_init((char *) filename, TRUE)) == NULL)
        return NULL;

    /* Enforce some limits.  */
    if (mdef->n_sen > BAD_SENID) {
        E_ERROR("Number of senones exceeds limit: %d > %d\n",
                mdef->n_sen, BAD_SENID);
        mdef_free(mdef);
        return NULL;
    }
    if (mdef->n_sseq > BAD_SSID) {
        E_ERROR("Number of senone sequences exceeds limit: %d > %d\n",
                mdef->n_sseq, BAD_SSID);
        mdef_free(mdef);
        return NULL;
    }
    /* We use uint8 for ciphones */
    if (mdef->n_ciphone > 255) {
        E_ERROR("Number of phones exceeds limit: %d > %d\n",
                mdef->n_ciphone, 255);
        mdef_free(mdef);
        return NULL;
    }

    bmdef = ckd_calloc(1, sizeof(*bmdef));
    bmdef->refcnt = 1;

    /* Easy stuff.  The mdef.c code has done the heavy lifting for us. */
    bmdef->n_ciphone = mdef->n_ciphone;
    bmdef->n_phone = mdef->n_phone;
    bmdef->n_emit_state = mdef->n_emit_state;
    bmdef->n_ci_sen = mdef->n_ci_sen;
    bmdef->n_sen = mdef->n_sen;
    bmdef->n_tmat = mdef->n_tmat;
    bmdef->n_sseq = mdef->n_sseq;
    bmdef->sseq = mdef->sseq;
    bmdef->cd2cisen = mdef->cd2cisen;
    bmdef->sen2cimap = mdef->sen2cimap;
    bmdef->n_ctx = 3;           /* Triphones only. */
    bmdef->sil = mdef->sil;
    mdef->sseq = NULL;          /* We are taking over this one. */
    mdef->cd2cisen = NULL;      /* And this one. */
    mdef->sen2cimap = NULL;     /* And this one. */

    /* Get the phone names.  If they are not sorted
     * ASCII-betically then we are in a world of hurt and
     * therefore will simply refuse to continue. */
    bmdef->ciname = ckd_calloc(bmdef->n_ciphone, sizeof(*bmdef->ciname));
    nchars = 0;
    for (i = 0; i < bmdef->n_ciphone; ++i)
        nchars += strlen(mdef->ciphone[i].name) + 1;
    bmdef->ciname[0] = ckd_calloc(nchars, 1);
    strcpy(bmdef->ciname[0], mdef->ciphone[0].name);
    for (i = 1; i < bmdef->n_ciphone; ++i) {
        bmdef->ciname[i] =
            bmdef->ciname[i - 1] + strlen(bmdef->ciname[i - 1]) + 1;
        strcpy(bmdef->ciname[i], mdef->ciphone[i].name);
        if (i > 0 && strcmp(bmdef->ciname[i - 1], bmdef->ciname[i]) > 0) {
            /* FIXME: there should be a solution to this, actually. */
            E_ERROR("Phone names are not in sorted order, sorry.");
            bin_mdef_free(bmdef);
            return NULL;
        }
    }

    /* Copy over phone information. */
    bmdef->phone = ckd_calloc(bmdef->n_phone, sizeof(*bmdef->phone));
    for (i = 0; i < mdef->n_phone; ++i) {
        bmdef->phone[i].ssid = mdef->phone[i].ssid;
        bmdef->phone[i].tmat = mdef->phone[i].tmat;
        if (i < bmdef->n_ciphone) {
            bmdef->phone[i].info.ci.filler = mdef->ciphone[i].filler;
        }
        else {
            bmdef->phone[i].info.cd.wpos = mdef->phone[i].wpos;
            bmdef->phone[i].info.cd.ctx[0] = mdef->phone[i].ci;
            bmdef->phone[i].info.cd.ctx[1] = mdef->phone[i].lc;
            bmdef->phone[i].info.cd.ctx[2] = mdef->phone[i].rc;
        }
    }

    /* Walk the wpos_ci_lclist once to find the total number of
     * nodes and the starting locations for each level. */
    nodes = lc_idx = ci_idx = rc_idx = 0;
    for (i = 0; i < N_WORD_POSN; ++i) {
        int j;
        for (j = 0; j < mdef->n_ciphone; ++j) {
            ph_lc_t *lc;

            for (lc = mdef->wpos_ci_lclist[i][j]; lc; lc = lc->next) {
                ph_rc_t *rc;
                for (rc = lc->rclist; rc; rc = rc->next) {
                    ++nodes;    /* RC node */
                }
                ++nodes;        /* LC node */
                ++rc_idx;       /* Start of RC nodes (after LC nodes) */
            }
            ++nodes;            /* CI node */
            ++lc_idx;           /* Start of LC nodes (after CI nodes) */
            ++rc_idx;           /* Start of RC nodes (after CI and LC nodes) */
        }
        ++nodes;                /* wpos node */
        ++ci_idx;               /* Start of CI nodes (after wpos nodes) */
        ++lc_idx;               /* Start of LC nodes (after CI nodes) */
        ++rc_idx;               /* STart of RC nodes (after wpos, CI, and LC nodes) */
    }
    E_INFO("Allocating %d * %d bytes (%d KiB) for CD tree\n",
           nodes, sizeof(*bmdef->cd_tree), 
           nodes * sizeof(*bmdef->cd_tree) / 1024);
    bmdef->n_cd_tree = nodes;
    bmdef->cd_tree = ckd_calloc(nodes, sizeof(*bmdef->cd_tree));
    for (i = 0; i < N_WORD_POSN; ++i) {
        int j;

        bmdef->cd_tree[i].ctx = i;
        bmdef->cd_tree[i].n_down = mdef->n_ciphone;
        bmdef->cd_tree[i].c.down = ci_idx;
#if 0
        E_INFO("%d => %c (%d@%d)\n",
               i, (WPOS_NAME)[i],
               bmdef->cd_tree[i].n_down, bmdef->cd_tree[i].c.down);
#endif

        /* Now we can build the rest of the tree. */
        for (j = 0; j < mdef->n_ciphone; ++j) {
            ph_lc_t *lc;

            bmdef->cd_tree[ci_idx].ctx = j;
            bmdef->cd_tree[ci_idx].c.down = lc_idx;
            for (lc = mdef->wpos_ci_lclist[i][j]; lc; lc = lc->next) {
                ph_rc_t *rc;

                bmdef->cd_tree[lc_idx].ctx = lc->lc;
                bmdef->cd_tree[lc_idx].c.down = rc_idx;
                for (rc = lc->rclist; rc; rc = rc->next) {
                    bmdef->cd_tree[rc_idx].ctx = rc->rc;
                    bmdef->cd_tree[rc_idx].n_down = 0;
                    bmdef->cd_tree[rc_idx].c.pid = rc->pid;
#if 0
                    E_INFO("%d => %s %s %s %c (%d@%d)\n",
                           rc_idx,
                           bmdef->ciname[j],
                           bmdef->ciname[lc->lc],
                           bmdef->ciname[rc->rc],
                           (WPOS_NAME)[i],
                           bmdef->cd_tree[rc_idx].n_down,
                           bmdef->cd_tree[rc_idx].c.down);
#endif

                    ++bmdef->cd_tree[lc_idx].n_down;
                    ++rc_idx;
                }
                /* If there are no triphones here,
                 * this is considered a leafnode, so
                 * set the pid to -1. */
                if (bmdef->cd_tree[lc_idx].n_down == 0)
                    bmdef->cd_tree[lc_idx].c.pid = -1;
#if 0
                E_INFO("%d => %s %s %c (%d@%d)\n",
                       lc_idx,
                       bmdef->ciname[j],
                       bmdef->ciname[lc->lc],
                       (WPOS_NAME)[i],
                       bmdef->cd_tree[lc_idx].n_down,
                       bmdef->cd_tree[lc_idx].c.down);
#endif

                ++bmdef->cd_tree[ci_idx].n_down;
                ++lc_idx;
            }

            /* As above, so below. */
            if (bmdef->cd_tree[ci_idx].n_down == 0)
                bmdef->cd_tree[ci_idx].c.pid = -1;
#if 0
            E_INFO("%d => %d=%s (%d@%d)\n",
                   ci_idx, j, bmdef->ciname[j],
                   bmdef->cd_tree[ci_idx].n_down,
                   bmdef->cd_tree[ci_idx].c.down);
#endif

            ++ci_idx;
        }
    }

    mdef_free(mdef);

    bmdef->alloc_mode = BIN_MDEF_FROM_TEXT;
    return bmdef;
}

bin_mdef_t *
bin_mdef_retain(bin_mdef_t *m)
{
    ++m->refcnt;
    return m;
}

int
bin_mdef_free(bin_mdef_t * m)
{
    if (m == NULL)
        return 0;
    if (--m->refcnt > 0)
        return m->refcnt;

    switch (m->alloc_mode) {
    case BIN_MDEF_FROM_TEXT:
        ckd_free(m->ciname[0]);
        ckd_free(m->sseq[0]);
        ckd_free(m->phone);
        ckd_free(m->cd_tree);
        break;
    case BIN_MDEF_IN_MEMORY:
        ckd_free(m->ciname[0]);
        break;
    case BIN_MDEF_ON_DISK:
        break;
    }
    if (m->filemap)
        mmio_file_unmap(m->filemap);
    ckd_free(m->cd2cisen);
    ckd_free(m->sen2cimap);
    ckd_free(m->ciname);
    ckd_free(m->sseq);
    ckd_free(m);
    return 0;
}

static const char format_desc[] =
    "BEGIN FILE FORMAT DESCRIPTION\n"
    "int32 n_ciphone;    /**< Number of base (CI) phones */\n"
    "int32 n_phone;	     /**< Number of base (CI) phones + (CD) triphones */\n"
    "int32 n_emit_state; /**< Number of emitting states per phone (0 if heterogeneous) */\n"
    "int32 n_ci_sen;     /**< Number of CI senones; these are the first */\n"
    "int32 n_sen;	     /**< Number of senones (CI+CD) */\n"
    "int32 n_tmat;	     /**< Number of transition matrices */\n"
    "int32 n_sseq;       /**< Number of unique senone sequences */\n"
    "int32 n_ctx;	     /**< Number of phones of context */\n"
    "int32 n_cd_tree;    /**< Number of nodes in CD tree structure */\n"
    "int32 sil;	     /**< CI phone ID for silence */\n"
    "char ciphones[][];  /**< CI phone strings (null-terminated) */\n"
    "char padding[];     /**< Padding to a 4-bytes boundary */\n"
    "struct { int16 ctx; int16 n_down; int32 pid/down } cd_tree[];\n"
    "struct { int32 ssid; int32 tmat; int8 attr[4] } phones[];\n"
    "int16 sseq[];       /**< Unique senone sequences */\n"
    "int8 sseq_len[];    /**< Number of states in each sseq (none if homogeneous) */\n"
    "END FILE FORMAT DESCRIPTION\n";

bin_mdef_t *
bin_mdef_read(cmd_ln_t *config, const char *filename)
{
    bin_mdef_t *m;
    FILE *fh;
    size_t tree_start;
    int32 val, i, swap, pos, end;
    int32 *sseq_size;
    int do_mmap;

    /* Try to read it as text first. */
    if ((m = bin_mdef_read_text(config, filename)) != NULL)
        return m;

    E_INFO("Reading binary model definition: %s\n", filename);
    if ((fh = fopen(filename, "rb")) == NULL)
        return NULL;

    if (fread(&val, 4, 1, fh) != 1) {
        fclose(fh);
        E_ERROR_SYSTEM("Failed to read byte-order marker from %s\n",
                       filename);
        return NULL;
    }
    swap = 0;
    if (val == BIN_MDEF_OTHER_ENDIAN) {
        swap = 1;
        E_INFO("Must byte-swap %s\n", filename);
    }
    if (fread(&val, 4, 1, fh) != 1) {
        fclose(fh);
        E_ERROR_SYSTEM("Failed to read version from %s\n", filename);
        return NULL;
    }
    if (swap)
        SWAP_INT32(&val);
    if (val > BIN_MDEF_FORMAT_VERSION) {
        E_ERROR("File format version %d for %s is newer than library\n",
                val, filename);
        fclose(fh);
        return NULL;
    }
    if (fread(&val, 4, 1, fh) != 1) {
        fclose(fh);
        E_ERROR_SYSTEM("Failed to read header length from %s\n", filename);
        return NULL;
    }
    if (swap)
        SWAP_INT32(&val);
    /* Skip format descriptor. */
    fseek(fh, val, SEEK_CUR);

    /* Finally allocate it. */
    m = ckd_calloc(1, sizeof(*m));
    m->refcnt = 1;

    /* Check these, to make gcc/glibc shut up. */
#define FREAD_SWAP32_CHK(dest)                                          \
    if (fread((dest), 4, 1, fh) != 1) {                                 \
        fclose(fh);                                                     \
        ckd_free(m);                                                    \
        E_ERROR_SYSTEM("Failed to read %s from %s\n", #dest, filename); \
        return NULL;                                                    \
    }                                                                   \
    if (swap) SWAP_INT32(dest);
    
    FREAD_SWAP32_CHK(&m->n_ciphone);
    FREAD_SWAP32_CHK(&m->n_phone);
    FREAD_SWAP32_CHK(&m->n_emit_state);
    FREAD_SWAP32_CHK(&m->n_ci_sen);
    FREAD_SWAP32_CHK(&m->n_sen);
    FREAD_SWAP32_CHK(&m->n_tmat);
    FREAD_SWAP32_CHK(&m->n_sseq);
    FREAD_SWAP32_CHK(&m->n_ctx);
    FREAD_SWAP32_CHK(&m->n_cd_tree);
    FREAD_SWAP32_CHK(&m->sil);

    /* CI names are first in the file. */
    m->ciname = ckd_calloc(m->n_ciphone, sizeof(*m->ciname));

    /* Decide whether to read in the whole file or mmap it. */
    do_mmap = config ? cmd_ln_boolean_r(config, "-mmap") : TRUE;
    if (swap) {
        E_WARN("-mmap specified, but mdef is other-endian.  Will not memory-map.\n");
        do_mmap = FALSE;
    } 
    /* Actually try to mmap it. */
    if (do_mmap) {
        m->filemap = mmio_file_read(filename);
        if (m->filemap == NULL)
            do_mmap = FALSE;
    }
    pos = ftell(fh);
    if (do_mmap) {
        /* Get the base pointer from the memory map. */
        m->ciname[0] = (char *)mmio_file_ptr(m->filemap) + pos;
        /* Success! */
        m->alloc_mode = BIN_MDEF_ON_DISK;
    }
    else {
        /* Read everything into memory. */
        m->alloc_mode = BIN_MDEF_IN_MEMORY;
        fseek(fh, 0, SEEK_END);
        end = ftell(fh);
        fseek(fh, pos, SEEK_SET);
        m->ciname[0] = ckd_malloc(end - pos);
        if (fread(m->ciname[0], 1, end - pos, fh) != end - pos)
            E_FATAL("Failed to read %d bytes of data from %s\n", end - pos, filename);
    }

    for (i = 1; i < m->n_ciphone; ++i)
        m->ciname[i] = m->ciname[i - 1] + strlen(m->ciname[i - 1]) + 1;

    /* Skip past the padding. */
    tree_start =
        m->ciname[i - 1] + strlen(m->ciname[i - 1]) + 1 - m->ciname[0];
    tree_start = (tree_start + 3) & ~3;
    m->cd_tree = (cd_tree_t *) (m->ciname[0] + tree_start);
    if (swap) {
        for (i = 0; i < m->n_cd_tree; ++i) {
            SWAP_INT16(&m->cd_tree[i].ctx);
            SWAP_INT16(&m->cd_tree[i].n_down);
            SWAP_INT32(&m->cd_tree[i].c.down);
        }
    }
    m->phone = (mdef_entry_t *) (m->cd_tree + m->n_cd_tree);
    if (swap) {
        for (i = 0; i < m->n_phone; ++i) {
            SWAP_INT32(&m->phone[i].ssid);
            SWAP_INT32(&m->phone[i].tmat);
        }
    }
    sseq_size = (int32 *) (m->phone + m->n_phone);
    if (swap)
        SWAP_INT32(sseq_size);
    m->sseq = ckd_calloc(m->n_sseq, sizeof(*m->sseq));
    m->sseq[0] = (uint16 *) (sseq_size + 1);
    if (swap) {
        for (i = 0; i < *sseq_size; ++i)
            SWAP_INT16(m->sseq[0] + i);
    }
    if (m->n_emit_state) {
        for (i = 1; i < m->n_sseq; ++i)
            m->sseq[i] = m->sseq[0] + i * m->n_emit_state;
    }
    else {
        m->sseq_len = (uint8 *) (m->sseq[0] + *sseq_size);
        for (i = 1; i < m->n_sseq; ++i)
            m->sseq[i] = m->sseq[i - 1] + m->sseq_len[i - 1];
    }

    /* Now build the CD-to-CI mappings using the senone sequences.
     * This is the only really accurate way to do it, though it is
     * still inaccurate in the case of heterogeneous topologies or
     * cross-state tying. */
    m->cd2cisen = (int16 *) ckd_malloc(m->n_sen * sizeof(*m->cd2cisen));
    m->sen2cimap = (int16 *) ckd_malloc(m->n_sen * sizeof(*m->sen2cimap));

    /* Default mappings (identity, none) */
    for (i = 0; i < m->n_ci_sen; ++i)
        m->cd2cisen[i] = i;
    for (; i < m->n_sen; ++i)
        m->cd2cisen[i] = -1;
    for (i = 0; i < m->n_sen; ++i)
        m->sen2cimap[i] = -1;
    for (i = 0; i < m->n_phone; ++i) {
        int32 j, ssid = m->phone[i].ssid;

        for (j = 0; j < bin_mdef_n_emit_state_phone(m, i); ++j) {
            int s = bin_mdef_sseq2sen(m, ssid, j);
            int ci = bin_mdef_pid2ci(m, i);
            /* Take the first one and warn if we have cross-state tying. */
            if (m->sen2cimap[s] == -1)
                m->sen2cimap[s] = ci;
            if (m->sen2cimap[s] != ci)
                E_WARN
                    ("Senone %d is shared between multiple base phones\n",
                     s);

            if (j > bin_mdef_n_emit_state_phone(m, ci))
                E_WARN("CD phone %d has fewer states than CI phone %d\n",
                       i, ci);
            else
                m->cd2cisen[s] =
                    bin_mdef_sseq2sen(m, m->phone[ci].ssid, j);
        }
    }

    /* Set the silence phone. */
    m->sil = bin_mdef_ciphone_id(m, S3_SILENCE_CIPHONE);

    E_INFO
        ("%d CI-phone, %d CD-phone, %d emitstate/phone, %d CI-sen, %d Sen, %d Sen-Seq\n",
         m->n_ciphone, m->n_phone - m->n_ciphone, m->n_emit_state,
         m->n_ci_sen, m->n_sen, m->n_sseq);
    fclose(fh);
    return m;
}

int
bin_mdef_write(bin_mdef_t * m, const char *filename)
{
    FILE *fh;
    int32 val, i;

    if ((fh = fopen(filename, "wb")) == NULL)
        return -1;

    /* Byteorder marker. */
    val = BIN_MDEF_NATIVE_ENDIAN;
    fwrite(&val, 1, 4, fh);
    /* Version. */
    val = BIN_MDEF_FORMAT_VERSION;
    fwrite(&val, 1, sizeof(val), fh);

    /* Round the format descriptor size up to a 4-byte boundary. */
    val = ((sizeof(format_desc) + 3) & ~3);
    fwrite(&val, 1, sizeof(val), fh);
    fwrite(format_desc, 1, sizeof(format_desc), fh);
    /* Pad it with zeros. */
    i = 0;
    fwrite(&i, 1, val - sizeof(format_desc), fh);

    /* Binary header things. */
    fwrite(&m->n_ciphone, 4, 1, fh);
    fwrite(&m->n_phone, 4, 1, fh);
    fwrite(&m->n_emit_state, 4, 1, fh);
    fwrite(&m->n_ci_sen, 4, 1, fh);
    fwrite(&m->n_sen, 4, 1, fh);
    fwrite(&m->n_tmat, 4, 1, fh);
    fwrite(&m->n_sseq, 4, 1, fh);
    fwrite(&m->n_ctx, 4, 1, fh);
    fwrite(&m->n_cd_tree, 4, 1, fh);
    /* Write this as a 32-bit value to preserve alignment for the
     * non-mmap case (we want things aligned both from the
     * beginning of the file and the beginning of the phone
     * strings). */
    val = m->sil;
    fwrite(&val, 4, 1, fh);

    /* Phone strings. */
    for (i = 0; i < m->n_ciphone; ++i)
        fwrite(m->ciname[i], 1, strlen(m->ciname[i]) + 1, fh);
    /* Pad with zeros. */
    val = (ftell(fh) + 3) & ~3;
    i = 0;
    fwrite(&i, 1, val - ftell(fh), fh);

    /* Write CD-tree */
    fwrite(m->cd_tree, sizeof(*m->cd_tree), m->n_cd_tree, fh);
    /* Write phones */
    fwrite(m->phone, sizeof(*m->phone), m->n_phone, fh);
    if (m->n_emit_state) {
        /* Write size of sseq */
        val = m->n_sseq * m->n_emit_state;
        fwrite(&val, 4, 1, fh);

        /* Write sseq */
        fwrite(m->sseq[0], sizeof(**m->sseq),
               m->n_sseq * m->n_emit_state, fh);
    }
    else {
        int32 n;

        /* Calcluate size of sseq */
        n = 0;
        for (i = 0; i < m->n_sseq; ++i)
            n += m->sseq_len[i];

        /* Write size of sseq */
        fwrite(&n, 4, 1, fh);

        /* Write sseq */
        fwrite(m->sseq[0], sizeof(**m->sseq), n, fh);

        /* Write sseq_len */
        fwrite(m->sseq_len, 1, m->n_sseq, fh);
    }
    fclose(fh);

    return 0;
}

int
bin_mdef_write_text(bin_mdef_t * m, const char *filename)
{
    FILE *fh;
    int p, i, n_total_state;

    if (strcmp(filename, "-") == 0)
        fh = stdout;
    else {
        if ((fh = fopen(filename, "w")) == NULL)
            return -1;
    }

    fprintf(fh, "0.3\n");
    fprintf(fh, "%d n_base\n", m->n_ciphone);
    fprintf(fh, "%d n_tri\n", m->n_phone - m->n_ciphone);
    if (m->n_emit_state)
        n_total_state = m->n_phone * (m->n_emit_state + 1);
    else {
        n_total_state = 0;
        for (i = 0; i < m->n_phone; ++i)
            n_total_state += m->sseq_len[m->phone[i].ssid] + 1;
    }
    fprintf(fh, "%d n_state_map\n", n_total_state);
    fprintf(fh, "%d n_tied_state\n", m->n_sen);
    fprintf(fh, "%d n_tied_ci_state\n", m->n_ci_sen);
    fprintf(fh, "%d n_tied_tmat\n", m->n_tmat);
    fprintf(fh, "#\n# Columns definitions\n");
    fprintf(fh, "#%4s %3s %3s %1s %6s %4s %s\n",
            "base", "lft", "rt", "p", "attrib", "tmat",
            "     ... state id's ...");

    for (p = 0; p < m->n_ciphone; p++) {
        int n_state;

        fprintf(fh, "%5s %3s %3s %1s", m->ciname[p], "-", "-", "-");

        if (bin_mdef_is_fillerphone(m, p))
            fprintf(fh, " %6s", "filler");
        else
            fprintf(fh, " %6s", "n/a");
        fprintf(fh, " %4d", m->phone[p].tmat);

        if (m->n_emit_state)
            n_state = m->n_emit_state;
        else
            n_state = m->sseq_len[m->phone[p].ssid];
        for (i = 0; i < n_state; i++) {
            fprintf(fh, " %6u", m->sseq[m->phone[p].ssid][i]);
        }
        fprintf(fh, " N\n");
    }


    for (; p < m->n_phone; p++) {
        int n_state;

        fprintf(fh, "%5s %3s %3s %c",
                m->ciname[m->phone[p].info.cd.ctx[0]],
                m->ciname[m->phone[p].info.cd.ctx[1]],
                m->ciname[m->phone[p].info.cd.ctx[2]],
                (WPOS_NAME)[m->phone[p].info.cd.wpos]);

        if (bin_mdef_is_fillerphone(m, p))
            fprintf(fh, " %6s", "filler");
        else
            fprintf(fh, " %6s", "n/a");
        fprintf(fh, " %4d", m->phone[p].tmat);


        if (m->n_emit_state)
            n_state = m->n_emit_state;
        else
            n_state = m->sseq_len[m->phone[p].ssid];
        for (i = 0; i < n_state; i++) {
            fprintf(fh, " %6u", m->sseq[m->phone[p].ssid][i]);
        }
        fprintf(fh, " N\n");
    }

    if (strcmp(filename, "-") != 0)
        fclose(fh);
    return 0;
}

int
bin_mdef_ciphone_id(bin_mdef_t * m, const char *ciphone)
{
    int low, mid, high;

    /* Exact binary search on m->ciphone */
    low = 0;
    high = m->n_ciphone;
    while (low < high) {
        int c;

        mid = (low + high) / 2;
        c = strcmp(ciphone, m->ciname[mid]);
        if (c == 0)
            return mid;
        else if (c > 0)
            low = mid + 1;
        else if (c < 0)
            high = mid;
    }
    return -1;
}

int
bin_mdef_ciphone_id_nocase(bin_mdef_t * m, const char *ciphone)
{
    int low, mid, high;

    /* Exact binary search on m->ciphone */
    low = 0;
    high = m->n_ciphone;
    while (low < high) {
        int c;

        mid = (low + high) / 2;
        c = strcmp_nocase(ciphone, m->ciname[mid]);
        if (c == 0)
            return mid;
        else if (c > 0)
            low = mid + 1;
        else if (c < 0)
            high = mid;
    }
    return -1;
}

const char *
bin_mdef_ciphone_str(bin_mdef_t * m, int32 ci)
{
    assert(m != NULL);
    assert(ci < m->n_ciphone);
    return m->ciname[ci];
}

int
bin_mdef_phone_id(bin_mdef_t * m, int32 ci, int32 lc, int32 rc, int32 wpos)
{
    cd_tree_t *cd_tree;
    int level, max;
    int16 ctx[4];

    assert(m);

    /* In the future, we might back off when context is not available,
     * but for now we'll just return the CI phone. */
    if (lc < 0 || rc < 0)
        return ci;

    assert((ci >= 0) && (ci < m->n_ciphone));
    assert((lc >= 0) && (lc < m->n_ciphone));
    assert((rc >= 0) && (rc < m->n_ciphone));
    assert((wpos >= 0) && (wpos < N_WORD_POSN));

    /* Create a context list, mapping fillers to silence. */
    ctx[0] = wpos;
    ctx[1] = ci;
    ctx[2] = (m->sil >= 0
              && m->phone[lc].info.ci.filler) ? m->sil : lc;
    ctx[3] = (m->sil >= 0
              && m->phone[rc].info.ci.filler) ? m->sil : rc;

    /* Walk down the cd_tree. */
    cd_tree = m->cd_tree;
    level = 0;                  /* What level we are on. */
    max = N_WORD_POSN;          /* Number of nodes on this level. */
    while (level < 4) {
        int i;

#if 0
        E_INFO("Looking for context %d=%s in %d at %d\n",
               ctx[level], m->ciname[ctx[level]],
               max, cd_tree - m->cd_tree);
#endif
        for (i = 0; i < max; ++i) {
#if 0
            E_INFO("Look at context %d=%s at %d\n",
                   cd_tree[i].ctx,
                   m->ciname[cd_tree[i].ctx], cd_tree + i - m->cd_tree);
#endif
            if (cd_tree[i].ctx == ctx[level])
                break;
        }
        if (i == max)
            return -1;
#if 0
        E_INFO("Found context %d=%s at %d, n_down=%d, down=%d\n",
               ctx[level], m->ciname[ctx[level]],
               cd_tree + i - m->cd_tree,
               cd_tree[i].n_down, cd_tree[i].c.down);
#endif
        /* Leaf node, stop here. */
        if (cd_tree[i].n_down == 0)
            return cd_tree[i].c.pid;

        /* Go down one level. */
        max = cd_tree[i].n_down;
        cd_tree = m->cd_tree + cd_tree[i].c.down;
        ++level;
    }
    /* We probably shouldn't get here. */
    return -1;
}

int
bin_mdef_phone_id_nearest(bin_mdef_t * m, int32 b, int32 l, int32 r, int32 pos)
{
    int p, tmppos;



    /* In the future, we might back off when context is not available,
     * but for now we'll just return the CI phone. */
    if (l < 0 || r < 0)
        return b;

    p = bin_mdef_phone_id(m, b, l, r, pos);
    if (p >= 0)
        return p;

    /* Exact triphone not found; backoff to other word positions */
    for (tmppos = 0; tmppos < N_WORD_POSN; tmppos++) {
        if (tmppos != pos) {
            p = bin_mdef_phone_id(m, b, l, r, tmppos);
            if (p >= 0)
                return p;
        }
    }

    /* Nothing yet; backoff to silence phone if non-silence filler context */
    /* In addition, backoff to silence phone on left/right if in beginning/end position */
    if (m->sil >= 0) {
        int newl = l, newr = r;
        if (m->phone[(int)l].info.ci.filler
            || pos == WORD_POSN_BEGIN || pos == WORD_POSN_SINGLE)
            newl = m->sil;
        if (m->phone[(int)r].info.ci.filler
            || pos == WORD_POSN_END || pos == WORD_POSN_SINGLE)
            newr = m->sil;
        if ((newl != l) || (newr != r)) {
            p = bin_mdef_phone_id(m, b, newl, newr, pos);
            if (p >= 0)
                return p;

            for (tmppos = 0; tmppos < N_WORD_POSN; tmppos++) {
                if (tmppos != pos) {
                    p = bin_mdef_phone_id(m, b, newl, newr, tmppos);
                    if (p >= 0)
                        return p;
                }
            }
        }
    }

    /* Nothing yet; backoff to base phone */
    return b;
}

int
bin_mdef_phone_str(bin_mdef_t * m, int pid, char *buf)
{
    char *wpos_name;

    assert(m);
    assert((pid >= 0) && (pid < m->n_phone));
    wpos_name = WPOS_NAME;

    buf[0] = '\0';
    if (pid < m->n_ciphone)
        sprintf(buf, "%s", bin_mdef_ciphone_str(m, pid));
    else {
        sprintf(buf, "%s %s %s %c",
                bin_mdef_ciphone_str(m, m->phone[pid].info.cd.ctx[0]),
                bin_mdef_ciphone_str(m, m->phone[pid].info.cd.ctx[1]),
                bin_mdef_ciphone_str(m, m->phone[pid].info.cd.ctx[2]),
                wpos_name[m->phone[pid].info.cd.wpos]);
    }
    return 0;
}
/* -*- c-file-style: "linux" -*- */
/* ====================================================================
 * Copyright (c) 2005 Carnegie Mellon University.  All rights 
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/**
 * @file bin_mdef.h
 * 
 * Binary format model definition files, with support for
 * heterogeneous topologies and variable-size N-phones
 *
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */
#ifndef __BIN_MDEF_H__
#define __BIN_MDEF_H__

#ifdef __cplusplus
extern "C" {
#if 0
}; /* Fool Emacs */
#endif
#endif /* __cplusplus */

/* SphinxBase headers. */
#include <sphinxbase/mmio.h>
#include <sphinxbase/cmd_ln.h>
#include <pocketsphinx_export.h>

#include "mdef.h"

#define BIN_MDEF_FORMAT_VERSION 1
/* Little-endian machines will write "BMDF" to disk, big-endian ones "FDMB". */
#define BIN_MDEF_NATIVE_ENDIAN 0x46444d42 /* 'BMDF' in little-endian order */
#define BIN_MDEF_OTHER_ENDIAN 0x424d4446  /* 'BMDF' in big-endian order */
#ifdef __GNUC__
#define ATTRIBUTE_PACKED __attribute__((packed))
#else
#define ATTRIBUTE_PACKED
#endif

/**
 * Phone entry (on-disk, 12 bytes)
 */
typedef struct mdef_entry_s mdef_entry_t;
struct mdef_entry_s {
	int32 ssid; /**< Senone sequence ID */
	int32 tmat; /**< Transition matrix ID */
	/* FIXME: is any of this actually necessary? */
	union {
		/**< CI phone information - attributes (just "filler" for now) */
		struct {
			uint8 filler;
			uint8 reserved[3];
		} ci;
		/**< CD phone information - context info. */
		struct {
			uint8 wpos;
			uint8 ctx[3]; /**< quintphones will require hacking */
		} cd;
	} info;
} ATTRIBUTE_PACKED;

/**
 * Invalid senone sequence ID (limited to 16 bits for PocketSphinx).
 */
#define BAD_SSID 0xffff
/**
 * Invalid senone ID (limited to 16 bits for PocketSphinx).
 */
#define BAD_SENID 0xffff

/**
 * Node in CD phone tree (on-disk, 8 bytes).
 */
typedef struct cd_tree_s cd_tree_t;
struct cd_tree_s {
	int16 ctx; /**< Context (word position or CI phone) */
	int16 n_down; /**< Number of children (0 for leafnode) */
	union {
		int32 pid; /**< Phone ID (leafnode) */
		int32 down; /**< Next level of the tree (offset from start of cd_trees) */ 
	} c;
};

/**
 * Model definition structure (in-memory).
 */
typedef struct bin_mdef_s bin_mdef_t;
struct bin_mdef_s {
	int refcnt;
	int32 n_ciphone;    /**< Number of base (CI) phones */
	int32 n_phone;	    /**< Number of base (CI) phones + (CD) triphones */
	int32 n_emit_state; /**< Number of emitting states per phone (0 for heterogeneous) */
	int32 n_ci_sen;	    /**< Number of CI senones; these are the first */
	int32 n_sen;	    /**< Number of senones (CI+CD) */
	int32 n_tmat;	    /**< Number of transition matrices */
	int32 n_sseq;       /**< Number of unique senone sequences */
	int32 n_ctx;	    /**< Number of phones of context */
	int32 n_cd_tree;    /**< Number of nodes in cd_tree (below) */
	int16 sil;	    /**< CI phone ID for silence */

	mmio_file_t *filemap;/**< File map for this file (if any) */
	char **ciname;       /**< CI phone names */
	cd_tree_t *cd_tree;  /**< Tree mapping CD phones to phone IDs */
	mdef_entry_t *phone; /**< All phone structures */
	uint16 **sseq;       /**< Unique senone sequences (2D array built at load time) */
	uint8 *sseq_len;     /**< Number of states in each sseq (NULL for homogeneous) */

	/* These two are not stored on disk, but are generated at load time. */
	int16 *cd2cisen;	/**< Parent CI-senone id for each senone */
	int16 *sen2cimap;	/**< Parent CI-phone for each senone (CI or CD) */

	/** Allocation mode for this object. */
	enum { BIN_MDEF_FROM_TEXT, BIN_MDEF_IN_MEMORY, BIN_MDEF_ON_DISK } alloc_mode;
};

#define bin_mdef_is_fillerphone(m,p)	(((p) < (m)->n_ciphone) \
		                         ? (m)->phone[p].info.ci.filler \
					 : (m)->phone[(m)->phone[p].info.cd.ctx[0]].info.ci.filler)
#define bin_mdef_is_ciphone(m,p)	((p) < (m)->n_ciphone)
#define bin_mdef_n_ciphone(m)		((m)->n_ciphone)
#define bin_mdef_n_phone(m)		((m)->n_phone)
#define bin_mdef_n_sseq(m)		((m)->n_sseq)
#define bin_mdef_n_emit_state(m)	((m)->n_emit_state)
#define bin_mdef_n_emit_state_phone(m,p) ((m)->n_emit_state ? (m)->n_emit_state \
					  : (m)->sseq_len[(m)->phone[p].ssid])
#define bin_mdef_n_sen(m)		((m)->n_sen)
#define bin_mdef_n_tmat(m)		((m)->n_tmat)
#define bin_mdef_pid2ssid(m,p)		((m)->phone[p].ssid)
#define bin_mdef_pid2tmatid(m,p)	((m)->phone[p].tmat)
#define bin_mdef_silphone(m)		((m)->sil)
#define bin_mdef_sen2cimap(m,s)		((m)->sen2cimap[s])
#define bin_mdef_sseq2sen(m,ss,pos)	((m)->sseq[ss][pos])
#define bin_mdef_pid2ci(m,p)		(((p) < (m)->n_ciphone) ? (p) \
                                         : (m)->phone[p].info.cd.ctx[0])

/**
 * Read a binary mdef from a file.
 */
POCKETSPHINX_EXPORT
bin_mdef_t *bin_mdef_read(cmd_ln_t *config, const char *filename);
/**
 * Read a text mdef from a file (creating an in-memory binary mdef).
 */
POCKETSPHINX_EXPORT
bin_mdef_t *bin_mdef_read_text(cmd_ln_t *config, const char *filename);
/**
 * Write a binary mdef to a file.
 */
POCKETSPHINX_EXPORT
int bin_mdef_write(bin_mdef_t *m, const char *filename);
/**
 * Write a binary mdef to a text file.
 */
POCKETSPHINX_EXPORT
int bin_mdef_write_text(bin_mdef_t *m, const char *filename);
/**
 * Retain a pointer to a bin_mdef_t.
 */
bin_mdef_t *bin_mdef_retain(bin_mdef_t *m);
/**
 * Release a pointer to a binary mdef.
 */
int bin_mdef_free(bin_mdef_t *m);

/**
 * Context-independent phone lookup.
 * @return phone id for ciphone.
 */
int bin_mdef_ciphone_id(bin_mdef_t *m,	       /**< In: Model structure being queried */
			const char *ciphone);  /**< In: ciphone for which id wanted */

/**
 * Case-insensitive context-independent phone lookup.
 * @return phone id for ciphone.
 */
int bin_mdef_ciphone_id_nocase(bin_mdef_t *m,	     /**< In: Model structure being queried */
			       const char *ciphone); /**< In: ciphone for which id wanted */

/* Return value: READ-ONLY ciphone string name for the given ciphone id */
const char *bin_mdef_ciphone_str(bin_mdef_t *m,	/**< In: Model structure being queried */
				 int32 ci);	/**< In: ciphone id for which name wanted */

/* Return value: phone id for the given constituents if found, else -1 */
int bin_mdef_phone_id(bin_mdef_t *m,	/**< In: Model structure being queried */
		      int32 b,		/**< In: base ciphone id */
		      int32 l,		/**< In: left context ciphone id */
		      int32 r,		/**< In: right context ciphone id */
		      int32 pos);	/**< In: Word position */

/* Look up a phone id, backing off to other word positions. */
int bin_mdef_phone_id_nearest(bin_mdef_t * m, int32 b,
			      int32 l, int32 r, int32 pos);

/**
 * Create a phone string for the given phone (base or triphone) id in the given buf.
 *
 * @return 0 if successful, -1 if error.
 */
int bin_mdef_phone_str(bin_mdef_t *m,	/**< In: Model structure being queried */
		       int pid,		/**< In: phone id being queried */
		       char *buf);	/**< Out: On return, buf has the string */

#ifdef __cplusplus
}; /* extern "C" */
#endif /* __cplusplus */

#endif /* __BIN_MDEF_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * blkarray_list.c -- block array-based list structure.
 * 
 * HISTORY
 * 
 * 18-Feb-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon
 * 		Started.
 */

/* System headers. */
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>

/* Local headers. */
#include "blkarray_list.h"


#define BLKARRAY_DEFAULT_MAXBLKS	16380
#define BLKARRAY_DEFAULT_BLKSIZE	16380


blkarray_list_t *
_blkarray_list_init(int32 maxblks, int32 blksize)
{
    blkarray_list_t *bl;

    if ((maxblks <= 0) || (blksize <= 0)) {
        E_ERROR("Cannot allocate %dx%d blkarray\n", maxblks, blksize);
        return NULL;
    }

    bl = (blkarray_list_t *) ckd_calloc(1, sizeof(blkarray_list_t));
    bl->ptr = (void ***) ckd_calloc(maxblks, sizeof(void **));
    bl->maxblks = maxblks;
    bl->blksize = blksize;
    bl->n_valid = 0;
    bl->cur_row = -1;           /* No row is allocated (dummy) */
    bl->cur_row_free = blksize; /* The dummy row is full */

    return bl;
}


blkarray_list_t *
blkarray_list_init(void)
{
    return _blkarray_list_init(BLKARRAY_DEFAULT_MAXBLKS,
                               BLKARRAY_DEFAULT_BLKSIZE);
}

void
blkarray_list_free(blkarray_list_t *bl)
{
    blkarray_list_reset(bl);
    ckd_free(bl->ptr);
    ckd_free(bl);
}


int32
blkarray_list_append(blkarray_list_t * bl, void *data)
{
    int32 id;

    assert(bl);

    if (bl->cur_row_free >= bl->blksize) {
        /* Previous row is filled; need to allocate a new row */
        bl->cur_row++;

        if (bl->cur_row >= bl->maxblks) {
            E_ERROR("Block array (%dx%d) exhausted\n",
                    bl->maxblks, bl->blksize);
            bl->cur_row--;
            return -1;
        }

        /* Allocate the new row */
        assert(bl->ptr[bl->cur_row] == NULL);
        bl->ptr[bl->cur_row] = (void **) ckd_malloc(bl->blksize *
                                                    sizeof(void *));

        bl->cur_row_free = 0;
    }

    bl->ptr[bl->cur_row][bl->cur_row_free] = data;
    (bl->cur_row_free)++;

    id = (bl->n_valid)++;
    assert(id >= 0);

    return id;
}


void
blkarray_list_reset(blkarray_list_t * bl)
{
    int32 i, j;

    /* Free all the allocated elements as well as the blocks */
    for (i = 0; i < bl->cur_row; i++) {
        for (j = 0; j < bl->blksize; j++)
            ckd_free(bl->ptr[i][j]);

        ckd_free(bl->ptr[i]);
        bl->ptr[i] = NULL;
    }
    if (i == bl->cur_row) {     /* NEED THIS! (in case cur_row < 0) */
        for (j = 0; j < bl->cur_row_free; j++)
            ckd_free(bl->ptr[i][j]);

        ckd_free(bl->ptr[i]);
        bl->ptr[i] = NULL;
    }

    bl->n_valid = 0;
    bl->cur_row = -1;
    bl->cur_row_free = bl->blksize;
}
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * blkarray_list.h -- array-based list structure, for memory and access
 * 	efficiency.
 * 
 * HISTORY
 * 
 * $Log: blkarray_list.h,v $
 * Revision 1.1.1.1  2006/05/23 18:45:02  dhuggins
 * re-importation
 *
 * Revision 1.2  2004/12/10 16:48:58  rkm
 * Added continuous density acoustic model handling
 *
 * Revision 1.1  2004/07/16 00:57:12  egouvea
 * Added Ravi's implementation of FSG support.
 *
 * Revision 1.2  2004/05/27 14:22:57  rkm
 * FSG cross-word triphones completed (but for single-phone words)
 *
 * Revision 1.1.1.1  2004/03/01 14:30:31  rkm
 *
 *
 * Revision 1.1  2004/02/26 01:14:48  rkm
 * *** empty log message ***
 *
 * 
 * 18-Feb-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon
 * 		Started.
 */


#ifndef __S2_BLKARRAY_LIST_H__
#define __S2_BLKARRAY_LIST_H__


#include <sphinxbase/prim_type.h>


/*
 * For maintaining a (conceptual) "list" of pointers to arbitrary data.
 * The application is responsible for knowing the true data type.
 * Use an array instead of a true list for efficiency (both memory and
 * speed).  But use a blocked (2-D) array to allow dynamic resizing at a
 * coarse grain.  An entire block is allocated or freed, as appropriate.
 */
typedef struct blkarray_list_s {
  void ***ptr;		/* ptr[][] is the user-supplied ptr */
  int32 maxblks;	/* size of ptr (#rows) */
  int32 blksize;	/* size of ptr[] (#cols, ie, size of each row) */
  int32 n_valid;	/* # entries actually stored in the list */
  int32 cur_row;	/* The current row being that has empty entry */
  int32 cur_row_free;	/* First entry valid within the current row */
} blkarray_list_t;

/* Access macros */
#define blkarray_list_ptr(l,r,c)	((l)->ptr[r][c])
#define blkarray_list_maxblks(l)	((l)->maxblks)
#define blkarray_list_blksize(l)	((l)->blksize)
#define blkarray_list_n_valid(l)	((l)->n_valid)
#define blkarray_list_cur_row(l)	((l)->cur_row)
#define blkarray_list_cur_row_free(l)	((l)->cur_row_free)


/*
 * Initialize and return a new blkarray_list containing an empty list
 * (i.e., 0 length).  Sized for the given values of maxblks and blksize.
 * NOTE: (maxblks * blksize) should not overflow int32, but this is not
 * checked.
 * Return the allocated entry if successful, NULL if any error.
 */
blkarray_list_t *_blkarray_list_init (int32 maxblks, int32 blksize);


/*
 * Like _blkarray_list_init() above, but for some default values of
 * maxblks and blksize.
 */
blkarray_list_t *blkarray_list_init ( void );

/**
 * Completely finalize a blkarray_list.
 */
void blkarray_list_free(blkarray_list_t *bl);


/*
 * Append the given new entry (data) to the end of the list.
 * Return the index of the entry if successful, -1 if any error.
 * The returned indices are guaranteed to be successive integers (i.e.,
 * 0, 1, 2...) for successive append operations, until the list is reset,
 * when they resume from 0.
 */
int32 blkarray_list_append (blkarray_list_t *, void *data);


/*
 * Free all the entries in the list (using ckd_free) and reset the
 * list length to 0.
 */
void blkarray_list_reset (blkarray_list_t *);


#endif
/* -*- c-basic-offset:4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

#include <string.h>

#include "dict2pid.h"
#include "hmm.h"


/**
 * @file dict2pid.c - dictionary word to senone sequence mappings
 */

void
compress_table(s3ssid_t * uncomp_tab, s3ssid_t * com_tab,
               s3cipid_t * ci_map, int32 n_ci)
{
    int32 found;
    int32 r;
    int32 tmp_r;

    for (r = 0; r < n_ci; r++) {
        com_tab[r] = BAD_S3SSID;
        ci_map[r] = BAD_S3CIPID;
    }
    /** Compress this map */
    for (r = 0; r < n_ci; r++) {

        found = 0;
        for (tmp_r = 0; tmp_r < r && com_tab[tmp_r] != BAD_S3SSID; tmp_r++) {   /* If it appears before, just filled in cimap; */
            if (uncomp_tab[r] == com_tab[tmp_r]) {
                found = 1;
                ci_map[r] = tmp_r;
                break;
            }
        }

        if (found == 0) {
            com_tab[tmp_r] = uncomp_tab[r];
            ci_map[r] = tmp_r;
        }
    }
}


static void
compress_right_context_tree(dict2pid_t * d2p,
                            s3ssid_t ***rdiph_rc)
{
    int32 n_ci;
    int32 b, l, r;
    s3ssid_t *rmap;
    s3ssid_t *tmpssid;
    s3cipid_t *tmpcimap;
    bin_mdef_t *mdef = d2p->mdef;
    size_t alloc;

    n_ci = mdef->n_ciphone;

    tmpssid = ckd_calloc(n_ci, sizeof(s3ssid_t));
    tmpcimap = ckd_calloc(n_ci, sizeof(s3cipid_t));

    d2p->rssid =
        (xwdssid_t **) ckd_calloc(mdef->n_ciphone, sizeof(xwdssid_t *));
    alloc = mdef->n_ciphone * sizeof(xwdssid_t *);

    for (b = 0; b < n_ci; b++) {
        d2p->rssid[b] =
            (xwdssid_t *) ckd_calloc(mdef->n_ciphone, sizeof(xwdssid_t));
        alloc += mdef->n_ciphone * sizeof(xwdssid_t);

        for (l = 0; l < n_ci; l++) {
            rmap = rdiph_rc[b][l];
            compress_table(rmap, tmpssid, tmpcimap, mdef->n_ciphone);

            for (r = 0; r < mdef->n_ciphone && tmpssid[r] != BAD_S3SSID;
                 r++);

            if (tmpssid[0] != BAD_S3SSID) {
                d2p->rssid[b][l].ssid = ckd_calloc(r, sizeof(s3ssid_t));
                memcpy(d2p->rssid[b][l].ssid, tmpssid,
                       r * sizeof(s3ssid_t));
                d2p->rssid[b][l].cimap =
                    ckd_calloc(mdef->n_ciphone, sizeof(s3cipid_t));
                memcpy(d2p->rssid[b][l].cimap, tmpcimap,
                       (mdef->n_ciphone) * sizeof(s3cipid_t));
                d2p->rssid[b][l].n_ssid = r;
            }
            else {
                d2p->rssid[b][l].ssid = NULL;
                d2p->rssid[b][l].cimap = NULL;
                d2p->rssid[b][l].n_ssid = 0;
            }
        }
    }

    E_INFO("Allocated %d bytes (%d KiB) for word-final triphones\n",
           (int)alloc, (int)alloc / 1024);
    ckd_free(tmpssid);
    ckd_free(tmpcimap);
}

static void
compress_left_right_context_tree(dict2pid_t * d2p)
{
    int32 n_ci;
    int32 b, l, r;
    s3ssid_t *rmap;
    s3ssid_t *tmpssid;
    s3cipid_t *tmpcimap;
    bin_mdef_t *mdef = d2p->mdef;
    size_t alloc;

    n_ci = mdef->n_ciphone;

    tmpssid = ckd_calloc(n_ci, sizeof(s3ssid_t));
    tmpcimap = ckd_calloc(n_ci, sizeof(s3cipid_t));

    assert(d2p->lrdiph_rc);

    d2p->lrssid =
        (xwdssid_t **) ckd_calloc(mdef->n_ciphone, sizeof(xwdssid_t *));
    alloc = mdef->n_ciphone * sizeof(xwdssid_t *);

    for (b = 0; b < n_ci; b++) {

        d2p->lrssid[b] =
            (xwdssid_t *) ckd_calloc(mdef->n_ciphone, sizeof(xwdssid_t));
        alloc += mdef->n_ciphone * sizeof(xwdssid_t);

        for (l = 0; l < n_ci; l++) {
            rmap = d2p->lrdiph_rc[b][l];

            compress_table(rmap, tmpssid, tmpcimap, mdef->n_ciphone);

            for (r = 0; r < mdef->n_ciphone && tmpssid[r] != BAD_S3SSID;
                 r++);

            if (tmpssid[0] != BAD_S3SSID) {
                d2p->lrssid[b][l].ssid = ckd_calloc(r, sizeof(s3ssid_t));
                memcpy(d2p->lrssid[b][l].ssid, tmpssid,
                       r * sizeof(s3ssid_t));
                d2p->lrssid[b][l].cimap =
                    ckd_calloc(mdef->n_ciphone, sizeof(s3cipid_t));
                memcpy(d2p->lrssid[b][l].cimap, tmpcimap,
                       (mdef->n_ciphone) * sizeof(s3cipid_t));
                d2p->lrssid[b][l].n_ssid = r;
            }
            else {
                d2p->lrssid[b][l].ssid = NULL;
                d2p->lrssid[b][l].cimap = NULL;
                d2p->lrssid[b][l].n_ssid = 0;
            }
        }
    }

    /* Try to compress lrdiph_rc into lrdiph_rc_compressed */
    ckd_free(tmpssid);
    ckd_free(tmpcimap);

    E_INFO("Allocated %d bytes (%d KiB) for single-phone word triphones\n",
           (int)alloc, (int)alloc / 1024);
}

/**
   ARCHAN, A duplicate of get_rc_npid in ctxt_table.h.  I doubt whether it is correct
   because the compressed map has not been checked. 
*/
int32
get_rc_nssid(dict2pid_t * d2p, s3wid_t w)
{
    int32 pronlen;
    s3cipid_t b, lc;
    dict_t *dict = d2p->dict;

    pronlen = dict->word[w].pronlen;
    b = dict->word[w].ciphone[pronlen - 1];

    if (pronlen == 1) {
        /* Is this true ?
           No known left context.  But all cimaps (for any l) are identical; pick one 
        */
        /*E_INFO("Single phone word\n"); */
        return (d2p->lrssid[b][0].n_ssid);
    }
    else {
        /*    E_INFO("Multiple phone word\n"); */
        lc = dict->word[w].ciphone[pronlen - 2];
        return (d2p->rssid[b][lc].n_ssid);
    }

}

s3cipid_t *
dict2pid_get_rcmap(dict2pid_t * d2p, s3wid_t w)
{
    int32 pronlen;
    s3cipid_t b, lc;
    dict_t *dict = d2p->dict;

    pronlen = dict->word[w].pronlen;
    b = dict->word[w].ciphone[pronlen - 1];

    if (pronlen == 1) {
        /* Is this true ?
           No known left context.  But all cimaps (for any l) are identical; pick one 
        */
        /*E_INFO("Single phone word\n"); */
        return (d2p->lrssid[b][0].cimap);
    }
    else {
        /*    E_INFO("Multiple phone word\n"); */
        lc = dict->word[w].ciphone[pronlen - 2];
        return (d2p->rssid[b][lc].cimap);
    }
}

static void
free_compress_map(xwdssid_t ** tree, int32 n_ci)
{
    int32 b, l;
    for (b = 0; b < n_ci; b++) {
        for (l = 0; l < n_ci; l++) {
            ckd_free(tree[b][l].ssid);
            ckd_free(tree[b][l].cimap);
        }
        ckd_free(tree[b]);
    }
    ckd_free(tree);
}

static void
populate_lrdiph(dict2pid_t *d2p, s3ssid_t ***rdiph_rc, s3cipid_t b)
{
    bin_mdef_t *mdef = d2p->mdef;
    s3cipid_t l, r;

    for (l = 0; l < bin_mdef_n_ciphone(mdef); l++) {
        for (r = 0; r < bin_mdef_n_ciphone(mdef); r++) {
            s3pid_t p;
            p = bin_mdef_phone_id_nearest(mdef, (s3cipid_t) b,
                                          (s3cipid_t) l,
                                          (s3cipid_t) r,
                                          WORD_POSN_SINGLE);
            d2p->lrdiph_rc[b][l][r]
                = bin_mdef_pid2ssid(mdef, p);
            if (r == bin_mdef_silphone(mdef))
                d2p->ldiph_lc[b][r][l]
                    = bin_mdef_pid2ssid(mdef, p);
            if (rdiph_rc && l == bin_mdef_silphone(mdef))
                rdiph_rc[b][l][r]
                    = bin_mdef_pid2ssid(mdef, p);
            assert(IS_S3SSID(bin_mdef_pid2ssid(mdef, p)));
            E_DEBUG(2,("%s(%s,%s) => %d / %d\n",
                       bin_mdef_ciphone_str(mdef, b),
                       bin_mdef_ciphone_str(mdef, l),
                       bin_mdef_ciphone_str(mdef, r),
                       p, bin_mdef_pid2ssid(mdef, p)));
        }
    }
}

int
dict2pid_add_word(dict2pid_t *d2p,
                  int32 wid)
{
    bin_mdef_t *mdef = d2p->mdef;
    dict_t *d = d2p->dict;

    if (dict_pronlen(d, wid) > 1) {
        s3cipid_t l;
        /* Make sure we have left and right context diphones for this
         * word. */
        if (d2p->ldiph_lc[dict_first_phone(d, wid)][dict_second_phone(d, wid)][0]
            == BAD_S3SSID) {
            E_DEBUG(2, ("Filling in left-context diphones for %s(?,%s)\n",
                   bin_mdef_ciphone_str(mdef, dict_first_phone(d, wid)),
                   bin_mdef_ciphone_str(mdef, dict_second_phone(d, wid))));
            for (l = 0; l < bin_mdef_n_ciphone(mdef); l++) {
                int p
                    = bin_mdef_phone_id_nearest(mdef,
                                                dict_first_phone(d, wid), l,
                                                dict_second_phone(d, wid),
                                                WORD_POSN_BEGIN);
                d2p->ldiph_lc[dict_first_phone(d, wid)][dict_second_phone(d, wid)][l]
                    = bin_mdef_pid2ssid(mdef, p);
            }
        }
        if (d2p->rssid[dict_last_phone(d, wid)][dict_second_last_phone(d, wid)].n_ssid
            == 0) {
            s3ssid_t *rmap;
            s3ssid_t *tmpssid;
            s3cipid_t *tmpcimap;
            s3cipid_t r;

            E_DEBUG(2, ("Filling in right-context diphones for %s(%s,?)\n",
                   bin_mdef_ciphone_str(mdef, dict_last_phone(d, wid)),
                   bin_mdef_ciphone_str(mdef, dict_second_last_phone(d, wid))));
            rmap = ckd_calloc(bin_mdef_n_ciphone(mdef), sizeof(*rmap));
            for (r = 0; r < bin_mdef_n_ciphone(mdef); r++) {
                int p
                    = bin_mdef_phone_id_nearest(mdef,
                                                dict_last_phone(d, wid),
                                                dict_second_last_phone(d, wid), r,
                                                WORD_POSN_END);
                rmap[r] = bin_mdef_pid2ssid(mdef, p);
            }
            tmpssid = ckd_calloc(bin_mdef_n_ciphone(mdef), sizeof(*tmpssid));
            tmpcimap = ckd_calloc(bin_mdef_n_ciphone(mdef), sizeof(*tmpcimap));
            compress_table(rmap, tmpssid, tmpcimap, bin_mdef_n_ciphone(mdef));
            for (r = 0; r < mdef->n_ciphone && tmpssid[r] != BAD_S3SSID; r++)
                ;
            d2p->rssid[dict_last_phone(d, wid)][dict_second_last_phone(d, wid)].ssid = tmpssid;
            d2p->rssid[dict_last_phone(d, wid)][dict_second_last_phone(d, wid)].cimap = tmpcimap;
            d2p->rssid[dict_last_phone(d, wid)][dict_second_last_phone(d, wid)].n_ssid = r;
            ckd_free(rmap);
        }
    }
    else {
        /* Make sure we have a left-right context triphone entry for
         * this word. */
        E_INFO("Filling in context triphones for %s(?,?)\n",
               bin_mdef_ciphone_str(mdef, dict_first_phone(d, wid)));
        if (d2p->lrdiph_rc[dict_first_phone(d, wid)][0][0] == BAD_S3SSID) {
            populate_lrdiph(d2p, NULL, dict_first_phone(d, wid));
        }
    }

    return 0;
}

s3ssid_t
dict2pid_internal(dict2pid_t *d2p,
                  int32 wid,
                  int pos)
{
    int b, l, r, p;
    dict_t *dict = d2p->dict;
    bin_mdef_t *mdef = d2p->mdef;

    if (pos == 0 || pos == dict_pronlen(dict, wid))
        return BAD_S3SSID;

    b = dict_pron(dict, wid, pos);
    l = dict_pron(dict, wid, pos - 1);
    r = dict_pron(dict, wid, pos + 1);
    p = bin_mdef_phone_id_nearest(mdef, (s3cipid_t) b,
                                  (s3cipid_t) l, (s3cipid_t) r,
                                  WORD_POSN_INTERNAL);
    return bin_mdef_pid2ssid(mdef, p);
}

dict2pid_t *
dict2pid_build(bin_mdef_t * mdef, dict_t * dict)
{
    dict2pid_t *dict2pid;
    s3ssid_t ***rdiph_rc;
    bitvec_t *ldiph, *rdiph, *single;
    int32 pronlen;
    int32 b, l, r, w, p;

    E_INFO("Building PID tables for dictionary\n");
    assert(mdef);
    assert(dict);

    dict2pid = (dict2pid_t *) ckd_calloc(1, sizeof(dict2pid_t));
    dict2pid->refcount = 1;
    dict2pid->mdef = bin_mdef_retain(mdef);
    dict2pid->dict = dict_retain(dict);
    E_INFO("Allocating %d^3 * %d bytes (%d KiB) for word-initial triphones\n",
           mdef->n_ciphone, sizeof(s3ssid_t),
           mdef->n_ciphone * mdef->n_ciphone * mdef->n_ciphone * sizeof(s3ssid_t) / 1024);
    dict2pid->ldiph_lc =
        (s3ssid_t ***) ckd_calloc_3d(mdef->n_ciphone, mdef->n_ciphone,
                                     mdef->n_ciphone, sizeof(s3ssid_t));
    /* Only used internally to generate rssid */
    rdiph_rc =
        (s3ssid_t ***) ckd_calloc_3d(mdef->n_ciphone, mdef->n_ciphone,
                                     mdef->n_ciphone, sizeof(s3ssid_t));

    dict2pid->lrdiph_rc = (s3ssid_t ***) ckd_calloc_3d(mdef->n_ciphone,
                                                       mdef->n_ciphone,
                                                       mdef->n_ciphone,
                                                       sizeof
                                                       (s3ssid_t));
    /* Actually could use memset for this, if BAD_S3SSID is guaranteed
     * to be 65535... */
    for (b = 0; b < mdef->n_ciphone; ++b) {
        for (r = 0; r < mdef->n_ciphone; ++r) {
            for (l = 0; l < mdef->n_ciphone; ++l) {
                dict2pid->ldiph_lc[b][r][l] = BAD_S3SSID;
                dict2pid->lrdiph_rc[b][l][r] = BAD_S3SSID;
                rdiph_rc[b][l][r] = BAD_S3SSID;
            }
        }
    }

    /* Track which diphones / ciphones have been seen. */
    ldiph = bitvec_alloc(mdef->n_ciphone * mdef->n_ciphone);
    rdiph = bitvec_alloc(mdef->n_ciphone * mdef->n_ciphone);
    single = bitvec_alloc(mdef->n_ciphone);

    for (w = 0; w < dict_size(dict2pid->dict); w++) {
        pronlen = dict_pronlen(dict, w);

        if (pronlen >= 2) {
            b = dict_first_phone(dict, w);
            r = dict_second_phone(dict, w);
            /* Populate ldiph_lc */
            if (bitvec_is_clear(ldiph, b * mdef->n_ciphone + r)) {
                /* Mark this diphone as done */
                bitvec_set(ldiph, b * mdef->n_ciphone + r);

                /* Record all possible ssids for b(?,r) */
                for (l = 0; l < bin_mdef_n_ciphone(mdef); l++) {
                    p = bin_mdef_phone_id_nearest(mdef, (s3cipid_t) b,
                                              (s3cipid_t) l, (s3cipid_t) r,
                                              WORD_POSN_BEGIN);
                    dict2pid->ldiph_lc[b][r][l] = bin_mdef_pid2ssid(mdef, p);
                }
            }


            /* Populate rdiph_rc */
            l = dict_second_last_phone(dict, w);
            b = dict_last_phone(dict, w);
            if (bitvec_is_clear(rdiph, b * mdef->n_ciphone + l)) {
                /* Mark this diphone as done */
                bitvec_set(rdiph, b * mdef->n_ciphone + l);

                for (r = 0; r < bin_mdef_n_ciphone(mdef); r++) {
                    p = bin_mdef_phone_id_nearest(mdef, (s3cipid_t) b,
                                              (s3cipid_t) l, (s3cipid_t) r,
                                              WORD_POSN_END);
                    rdiph_rc[b][l][r] = bin_mdef_pid2ssid(mdef, p);
                }
            }
        }
        else if (pronlen == 1) {
            b = dict_pron(dict, w, 0);
            E_DEBUG(1,("Building tables for single phone word %s phone %d = %s\n",
                       dict_wordstr(dict, w), b, bin_mdef_ciphone_str(mdef, b)));
            /* Populate lrdiph_rc (and also ldiph_lc, rdiph_rc if needed) */
            if (bitvec_is_clear(single, b)) {
                populate_lrdiph(dict2pid, rdiph_rc, b);
                bitvec_set(single, b);
            }
        }
    }

    bitvec_free(ldiph);
    bitvec_free(rdiph);
    bitvec_free(single);

    /* Try to compress rdiph_rc into rdiph_rc_compressed */
    compress_right_context_tree(dict2pid, rdiph_rc);
    compress_left_right_context_tree(dict2pid);

    ckd_free_3d(rdiph_rc);

    dict2pid_report(dict2pid);
    return dict2pid;
}

dict2pid_t *
dict2pid_retain(dict2pid_t *d2p)
{
    ++d2p->refcount;
    return d2p;
}

int
dict2pid_free(dict2pid_t * d2p)
{
    if (d2p == NULL)
        return 0;
    if (--d2p->refcount > 0)
        return d2p->refcount;

    if (d2p->ldiph_lc)
        ckd_free_3d((void ***) d2p->ldiph_lc);

    if (d2p->lrdiph_rc)
        ckd_free_3d((void ***) d2p->lrdiph_rc);

    if (d2p->rssid)
        free_compress_map(d2p->rssid, bin_mdef_n_ciphone(d2p->mdef));

    if (d2p->lrssid)
        free_compress_map(d2p->lrssid, bin_mdef_n_ciphone(d2p->mdef));

    bin_mdef_free(d2p->mdef);
    dict_free(d2p->dict);
    ckd_free(d2p);
    return 0;
}

void
dict2pid_report(dict2pid_t * d2p)
{
}

void
dict2pid_dump(FILE * fp, dict2pid_t * d2p)
{
    int32 w, p, pronlen;
    int32 i, j, b, l, r;
    bin_mdef_t *mdef = d2p->mdef;
    dict_t *dict = d2p->dict;

    fprintf(fp, "# INTERNAL (wd comssid ssid ssid ... ssid comssid)\n");
    for (w = 0; w < dict_size(dict); w++) {
        fprintf(fp, "%30s ", dict_wordstr(dict, w));

        pronlen = dict_pronlen(dict, w);
        for (p = 0; p < pronlen; p++)
            fprintf(fp, " %5d", dict2pid_internal(d2p, w, p));
        fprintf(fp, "\n");
    }
    fprintf(fp, "#\n");

    fprintf(fp, "# LDIPH_LC (b r l ssid)\n");
    for (b = 0; b < bin_mdef_n_ciphone(mdef); b++) {
        for (r = 0; r < bin_mdef_n_ciphone(mdef); r++) {
            for (l = 0; l < bin_mdef_n_ciphone(mdef); l++) {
                if (IS_S3SSID(d2p->ldiph_lc[b][r][l]))
                    fprintf(fp, "%6s %6s %6s %5d\n", bin_mdef_ciphone_str(mdef, (s3cipid_t) b), bin_mdef_ciphone_str(mdef, (s3cipid_t) r), bin_mdef_ciphone_str(mdef, (s3cipid_t) l), d2p->ldiph_lc[b][r][l]);      /* RAH, ldiph_lc is returning an int32, %d expects an int16 */
            }
        }
    }
    fprintf(fp, "#\n");

    fprintf(fp, "# SSEQ %d (senid senid ...)\n", mdef->n_sseq);
    for (i = 0; i < mdef->n_sseq; i++) {
        fprintf(fp, "%5d ", i);
        for (j = 0; j < bin_mdef_n_emit_state(mdef); j++)
            fprintf(fp, " %5d", mdef->sseq[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "#\n");
    fprintf(fp, "# END\n");

    fflush(fp);
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * dict2pid.h -- Triphones for dictionary
 * 
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1999 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * $Log$
 * Revision 1.1  2006/04/05  20:27:30  dhdfu
 * A Great Reorganzation of header files and executables
 * 
 * Revision 1.9  2006/02/22 21:05:16  arthchan2003
 * Merged from branch SPHINX3_5_2_RCI_IRII_BRANCH:
 *
 * 1, Added logic to handle bothe composite and non composite left
 * triphone.  Composite left triphone's logic (the original one) is
 * tested thoroughly. The non-composite triphone (or full triphone) is
 * found to have bugs.  The latter is fended off from the users in the
 * library level.
 *
 * 2, Fixed dox-doc.
 *
 * Revision 1.8.4.5  2005/11/17 06:13:49  arthchan2003
 * Use compressed right context in expansion in triphones.
 *
 * Revision 1.8.4.4  2005/10/17 04:48:45  arthchan2003
 * Free resource correctly in dict2pid.
 *
 * Revision 1.8.4.3  2005/10/07 19:03:38  arthchan2003
 * Added xwdssid_t structure.  Also added compression routines.
 *
 * Revision 1.8.4.2  2005/09/25 19:13:31  arthchan2003
 * Added optional full triphone expansion support when building context phone mapping.
 *
 * Revision 1.8.4.1  2005/07/17 05:20:30  arthchan2003
 * Fixed dox-doc.
 *
 * Revision 1.8  2005/06/21 21:03:49  arthchan2003
 * 1, Introduced a reporting routine. 2, Fixed doyxgen documentation, 3, Added  keyword.
 *
 * Revision 1.5  2005/06/13 04:02:57  archan
 * Fixed most doxygen-style documentation under libs3decoder.
 *
 * Revision 1.4  2005/04/21 23:50:26  archan
 * Some more refactoring on the how reporting of structures inside kbcore_t is done, it is now 50% nice. Also added class-based LM test case into test-decode.sh.in.  At this moment, everything in search mode 5 is already done.  It is time to test the idea whether the search can really be used.
 *
 * Revision 1.3  2005/03/30 01:22:46  archan
 * Fixed mistakes in last updates. Add
 *
 * 
 * 14-Sep-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Added dict2pid_comsseq2sen_active().
 * 
 * 04-May-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Started.
 */


#ifndef _S3_DICT2PID_H_
#define _S3_DICT2PID_H_

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/logmath.h>
#include <sphinxbase/bitvec.h>

/* Local headers. */
#include "s3types.h"
#include "bin_mdef.h"
#include "dict.h"

/** \file dict2pid.h
 * \brief Building triphones for a dictionary. 
 *
 * This is one of the more complicated parts of a cross-word
 * triphone model decoder.  The first and last phones of each word
 * get their left and right contexts, respectively, from other
 * words.  For single-phone words, both its contexts are from other
 * words, simultaneously.  As these words are not known beforehand,
 * life gets complicated.
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

/**
 * \struct xwdssid_t
 * \brief cross word triphone model structure 
 */

typedef struct {
    s3ssid_t  *ssid;	/**< Senone Sequence ID list for all context ciphones */
    s3cipid_t *cimap;	/**< Index into ssid[] above for each ci phone */
    int32     n_ssid;	/**< #Unique ssid in above, compressed ssid list */
} xwdssid_t;

/**
   \struct dict2pid_t
   \brief Building composite triphone (as well as word internal triphones) with the dictionary. 
*/

typedef struct {
    int refcount;

    bin_mdef_t *mdef;           /**< Model definition, used to generate
                                   internal ssids on the fly. */
    dict_t *dict;               /**< Dictionary this table refers to. */

    /*Notice the order of the arguments */
    /* FIXME: This is crying out for compression - in Mandarin we have
     * 180 context independent phones, which makes this an 11MB
     * array. */
    s3ssid_t ***ldiph_lc;	/**< For multi-phone words, [base][rc][lc] -> ssid; filled out for
				   word-initial base x rc combinations in current vocabulary */


    xwdssid_t **rssid;          /**< Right context state sequence id table 
                                   First dimension: base phone,
                                   Second dimension: left context. 
                                */


    s3ssid_t ***lrdiph_rc;      /**< For single-phone words, [base][lc][rc] -> ssid; filled out for
                                   single-phone base x lc combinations in current vocabulary */

    xwdssid_t **lrssid;          /**< Left-Right context state sequence id table 
                                    First dimension: base phone,
                                    Second dimension: left context. 
                                 */
} dict2pid_t;

/** Access macros; not designed for arbitrary use */
#define dict2pid_rssid(d,ci,lc)  (&(d)->rssid[ci][lc])
#define dict2pid_ldiph_lc(d,b,r,l) ((d)->ldiph_lc[b][r][l])
#define dict2pid_lrdiph_rc(d,b,l,r) ((d)->lrdiph_rc[b][l][r])

/**
 * Build the dict2pid structure for the given model/dictionary
 */
dict2pid_t *dict2pid_build(bin_mdef_t *mdef,   /**< A  model definition*/
                           dict_t *dict        /**< An initialized dictionary */
    );

/**
 * Retain a pointer to dict2pid
 */
dict2pid_t *dict2pid_retain(dict2pid_t *d2p);  

/**
 * Free the memory dict2pid structure
 */
int dict2pid_free(dict2pid_t *d2p /**< In: the d2p */
    );

/**
 * Return the senone sequence ID for the given word position.
 */
s3ssid_t dict2pid_internal(dict2pid_t *d2p,
                           int32 wid,
                           int pos);

/**
 * Add a word to the dict2pid structure (after adding it to dict).
 */
int dict2pid_add_word(dict2pid_t *d2p,
                      int32 wid);

/**
 * For debugging
 */
void dict2pid_dump(FILE *fp,        /**< In: a file pointer */
                   dict2pid_t *d2p /**< In: a dict2pid_t structure */
    );

/** Report a dict2pid data structure */
void dict2pid_report(dict2pid_t *d2p /**< In: a dict2pid_t structure */
    );

/**
 * Get number of rc 
 */
int32 get_rc_nssid(dict2pid_t *d2p,  /**< In: a dict2pid */
		   s3wid_t w         /**< In: a wid */
    );

/**
 * Get RC map 
 */
s3cipid_t* dict2pid_get_rcmap(dict2pid_t *d2p,  /**< In: a dict2pid */
			      s3wid_t w        /**< In: a wid */
    );

#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif


#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/* System headers. */
#include <string.h>

/* SphinxBase headers. */
#include <sphinxbase/pio.h>
#include <sphinxbase/strfuncs.h>

/* Local headers. */
#include "dict.h"


#define DELIM	" \t\n"         /* Set of field separator characters */
#define DEFAULT_NUM_PHONE	(MAX_S3CIPID+1)

#if WIN32
#define snprintf sprintf_s
#endif 

extern const char *const cmu6_lts_phone_table[];

static s3cipid_t
dict_ciphone_id(dict_t * d, const char *str)
{
    if (d->nocase)
        return bin_mdef_ciphone_id_nocase(d->mdef, str);
    else
        return bin_mdef_ciphone_id(d->mdef, str);
}


const char *
dict_ciphone_str(dict_t * d, s3wid_t wid, int32 pos)
{
    assert(d != NULL);
    assert((wid >= 0) && (wid < d->n_word));
    assert((pos >= 0) && (pos < d->word[wid].pronlen));

    return bin_mdef_ciphone_str(d->mdef, d->word[wid].ciphone[pos]);
}


s3wid_t
dict_add_word(dict_t * d, char const *word, s3cipid_t const * p, int32 np)
{
    int32 len;
    dictword_t *wordp;
    s3wid_t newwid;
    char *wword;

    if (d->n_word >= d->max_words) {
        E_INFO("Reallocating to %d KiB for word entries\n",
               (d->max_words + S3DICT_INC_SZ) * sizeof(dictword_t) / 1024);
        d->word =
            (dictword_t *) ckd_realloc(d->word,
                                       (d->max_words +
                                        S3DICT_INC_SZ) * sizeof(dictword_t));
        d->max_words = d->max_words + S3DICT_INC_SZ;
    }

    wordp = d->word + d->n_word;
    wordp->word = (char *) ckd_salloc(word);    /* Freed in dict_free */

    /* Associate word string with d->n_word in hash table */
    if (hash_table_enter_int32(d->ht, wordp->word, d->n_word) != d->n_word) {
        ckd_free(wordp->word);
        wordp->word = NULL;
        return BAD_S3WID;
    }

    /* Fill in word entry, and set defaults */
    if (p && (np > 0)) {
        wordp->ciphone = (s3cipid_t *) ckd_malloc(np * sizeof(s3cipid_t));      /* Freed in dict_free */
        memcpy(wordp->ciphone, p, np * sizeof(s3cipid_t));
        wordp->pronlen = np;
    }
    else {
        wordp->ciphone = NULL;
        wordp->pronlen = 0;
    }
    wordp->alt = BAD_S3WID;
    wordp->basewid = d->n_word;

    /* Determine base/alt wids */
    wword = ckd_salloc(word);
    if ((len = dict_word2basestr(wword)) > 0) {
	int32 w;

        /* Truncated to a baseword string; find its ID */
        if (hash_table_lookup_int32(d->ht, wword, &w) < 0) {
            E_ERROR("Missing base word for: %s\n", word);
            ckd_free(wword);
            ckd_free(wordp->word);
            wordp->word = NULL;
            return BAD_S3WID;
        }

        /* Link into alt list */
        wordp->basewid = w;
        wordp->alt = d->word[w].alt;
        d->word[w].alt = d->n_word;
    }
    ckd_free(wword);

    newwid = d->n_word++;

    return newwid;
}


static int32
dict_read(FILE * fp, dict_t * d)
{
    lineiter_t *li;
    char **wptr;
    s3cipid_t *p;
    int32 lineno, nwd;
    s3wid_t w;
    int32 i, maxwd;
    size_t stralloc, phnalloc;

    maxwd = 512;
    p = (s3cipid_t *) ckd_calloc(maxwd + 4, sizeof(*p));
    wptr = (char **) ckd_calloc(maxwd, sizeof(char *)); /* Freed below */

    lineno = 0;
    stralloc = phnalloc = 0;
    for (li = lineiter_start(fp); li; li = lineiter_next(li)) {
        lineno++;
        if (0 == strncmp(li->buf, "##", 2)
            || 0 == strncmp(li->buf, ";;", 2))
            continue;

        if ((nwd = str2words(li->buf, wptr, maxwd)) < 0) {
            /* Increase size of p, wptr. */
            nwd = str2words(li->buf, NULL, 0);
            assert(nwd > maxwd); /* why else would it fail? */
            maxwd = nwd;
            p = (s3cipid_t *) ckd_realloc(p, (maxwd + 4) * sizeof(*p));
            wptr = (char **) ckd_realloc(wptr, maxwd * sizeof(*wptr));
        }

        if (nwd == 0)           /* Empty line */
            continue;
        /* wptr[0] is the word-string and wptr[1..nwd-1] the pronunciation sequence */
        if (nwd == 1) {
            E_ERROR("Line %d: No pronunciation for word '%s'; ignored\n",
                    lineno, wptr[0]);
            continue;
        }


        /* Convert pronunciation string to CI-phone-ids */
        for (i = 1; i < nwd; i++) {
            p[i - 1] = dict_ciphone_id(d, wptr[i]);
            if (NOT_S3CIPID(p[i - 1])) {
                E_ERROR("Line %d: Phone '%s' is mising in the acoustic model; word '%s' ignored\n",
                        lineno, wptr[i], wptr[0]);
                break;
            }
        }

        if (i == nwd) {         /* All CI-phones successfully converted to IDs */
            w = dict_add_word(d, wptr[0], p, nwd - 1);
            if (NOT_S3WID(w))
                E_ERROR
                    ("Line %d: Failed to add the word '%s' (duplicate?); ignored\n",
                     lineno, wptr[0]);
            else {
                stralloc += strlen(d->word[w].word);
                phnalloc += d->word[w].pronlen * sizeof(s3cipid_t);
            }
        }
    }
    E_INFO("Allocated %d KiB for strings, %d KiB for phones\n",
           (int)stralloc / 1024, (int)phnalloc / 1024);
    ckd_free(p);
    ckd_free(wptr);

    return 0;
}

int
dict_write(dict_t *dict, char const *filename, char const *format)
{
    FILE *fh;
    int i;

    if ((fh = fopen(filename, "w")) == NULL) {
        E_ERROR_SYSTEM("Failed to open '%s'", filename);
        return -1;
    }
    for (i = 0; i < dict->n_word; ++i) {
        char *phones;
        int j, phlen;
        if (!dict_real_word(dict, i))
            continue;
        for (phlen = j = 0; j < dict_pronlen(dict, i); ++j)
            phlen += strlen(dict_ciphone_str(dict, i, j)) + 1;
        phones = ckd_calloc(1, phlen);
        for (j = 0; j < dict_pronlen(dict, i); ++j) {
            strcat(phones, dict_ciphone_str(dict, i, j));
            if (j != dict_pronlen(dict, i) - 1)
                strcat(phones, " ");
        }
        fprintf(fh, "%-30s %s\n", dict_wordstr(dict, i), phones);
        ckd_free(phones);
    }
    fclose(fh);
    return 0;
}


dict_t *
dict_init(cmd_ln_t *config, bin_mdef_t * mdef)
{
    FILE *fp, *fp2;
    int32 n;
    lineiter_t *li;
    dict_t *d;
    s3cipid_t sil;
    char const *dictfile = NULL, *fillerfile = NULL;

    if (config) {
        dictfile = cmd_ln_str_r(config, "-dict");
        fillerfile = cmd_ln_str_r(config, "-fdict");
    }

    /*
     * First obtain #words in dictionary (for hash table allocation).
     * Reason: The PC NT system doesn't like to grow memory gradually.  Better to allocate
     * all the required memory in one go.
     */
    fp = NULL;
    n = 0;
    if (dictfile) {
        if ((fp = fopen(dictfile, "r")) == NULL) {
            E_ERROR_SYSTEM("Failed to open dictionary file '%s' for reading", dictfile);
    	    return NULL;
        }
        for (li = lineiter_start(fp); li; li = lineiter_next(li)) {
	    if (0 != strncmp(li->buf, "##", 2)
    	        && 0 != strncmp(li->buf, ";;", 2))
                n++;
        }
        rewind(fp);
    }

    fp2 = NULL;
    if (fillerfile) {
        if ((fp2 = fopen(fillerfile, "r")) == NULL) {
            E_ERROR_SYSTEM("Failed to open filler dictionary file '%s' for reading", fillerfile);
            fclose(fp);
            return NULL;
	}
        for (li = lineiter_start(fp2); li; li = lineiter_next(li)) {
	    if (0 != strncmp(li->buf, "##", 2)
    	        && 0 != strncmp(li->buf, ";;", 2))
                n++;
        }
        rewind(fp2);
    }

    /*
     * Allocate dict entries.  HACK!!  Allow some extra entries for words not in file.
     * Also check for type size restrictions.
     */
    d = (dict_t *) ckd_calloc(1, sizeof(dict_t));       /* freed in dict_free() */
    d->refcnt = 1;
    d->max_words =
        (n + S3DICT_INC_SZ < MAX_S3WID) ? n + S3DICT_INC_SZ : MAX_S3WID;
    if (n >= MAX_S3WID) {
        E_ERROR("Number of words in dictionaries (%d) exceeds limit (%d)\n", n,
                MAX_S3WID);
        fclose(fp);
        fclose(fp2);
        ckd_free(d);
        return NULL;
    }

    E_INFO("Allocating %d * %d bytes (%d KiB) for word entries\n",
           d->max_words, sizeof(dictword_t),
           d->max_words * sizeof(dictword_t) / 1024);
    d->word = (dictword_t *) ckd_calloc(d->max_words, sizeof(dictword_t));      /* freed in dict_free() */
    d->n_word = 0;
    if (mdef)
        d->mdef = bin_mdef_retain(mdef);

    /* Create new hash table for word strings; case-insensitive word strings */
    if (config && cmd_ln_exists_r(config, "-dictcase"))
        d->nocase = cmd_ln_boolean_r(config, "-dictcase");
    d->ht = hash_table_new(d->max_words, d->nocase);

    /* Digest main dictionary file */
    if (fp) {
        E_INFO("Reading main dictionary: %s\n", dictfile);
        dict_read(fp, d);
        fclose(fp);
        E_INFO("%d words read\n", d->n_word);
    }

    /* Now the filler dictionary file, if it exists */
    d->filler_start = d->n_word;
    if (fillerfile) {
        E_INFO("Reading filler dictionary: %s\n", fillerfile);
        dict_read(fp2, d);
        fclose(fp2);
        E_INFO("%d words read\n", d->n_word - d->filler_start);
    }
    if (mdef)
        sil = bin_mdef_silphone(mdef);
    else
        sil = 0;
    if (dict_wordid(d, S3_START_WORD) == BAD_S3WID) {
        dict_add_word(d, S3_START_WORD, &sil, 1);
    }
    if (dict_wordid(d, S3_FINISH_WORD) == BAD_S3WID) {
        dict_add_word(d, S3_FINISH_WORD, &sil, 1);
    }
    if (dict_wordid(d, S3_SILENCE_WORD) == BAD_S3WID) {
        dict_add_word(d, S3_SILENCE_WORD, &sil, 1);
    }

    d->filler_end = d->n_word - 1;

    /* Initialize distinguished word-ids */
    d->startwid = dict_wordid(d, S3_START_WORD);
    d->finishwid = dict_wordid(d, S3_FINISH_WORD);
    d->silwid = dict_wordid(d, S3_SILENCE_WORD);

    if ((d->filler_start > d->filler_end)
        || (!dict_filler_word(d, d->silwid))) {
        E_ERROR("Word '%s' must occur (only) in filler dictionary\n",
                S3_SILENCE_WORD);
        dict_free(d);
        return NULL;
    }

    /* No check that alternative pronunciations for filler words are in filler range!! */

    return d;
}


s3wid_t
dict_wordid(dict_t *d, const char *word)
{
    int32 w;

    assert(d);
    assert(word);

    if (hash_table_lookup_int32(d->ht, word, &w) < 0)
        return (BAD_S3WID);
    return w;
}


int
dict_filler_word(dict_t *d, s3wid_t w)
{
    assert(d);
    assert((w >= 0) && (w < d->n_word));

    w = dict_basewid(d, w);
    if ((w == d->startwid) || (w == d->finishwid))
        return 0;
    if ((w >= d->filler_start) && (w <= d->filler_end))
        return 1;
    return 0;
}

int
dict_real_word(dict_t *d, s3wid_t w)
{
    assert(d);
    assert((w >= 0) && (w < d->n_word));

    w = dict_basewid(d, w);
    if ((w == d->startwid) || (w == d->finishwid))
        return 0;
    if ((w >= d->filler_start) && (w <= d->filler_end))
        return 0;
    return 1;
}


int32
dict_word2basestr(char *word)
{
    int32 i, len;

    len = strlen(word);
    if (word[len - 1] == ')') {
        for (i = len - 2; (i > 0) && (word[i] != '('); --i);

        if (i > 0) {
            /* The word is of the form <baseword>(...); strip from left-paren */
            word[i] = '\0';
            return i;
        }
    }

    return -1;
}

dict_t *
dict_retain(dict_t *d)
{
    ++d->refcnt;
    return d;
}

int
dict_free(dict_t * d)
{
    int i;
    dictword_t *word;

    if (d == NULL)
        return 0;
    if (--d->refcnt > 0)
        return d->refcnt;

    /* First Step, free all memory allocated for each word */
    for (i = 0; i < d->n_word; i++) {
        word = (dictword_t *) & (d->word[i]);
        if (word->word)
            ckd_free((void *) word->word);
        if (word->ciphone)
            ckd_free((void *) word->ciphone);
    }

    if (d->word)
        ckd_free((void *) d->word);
    if (d->ht)
        hash_table_free(d->ht);
    if (d->mdef)
        bin_mdef_free(d->mdef);
    ckd_free((void *) d);

    return 0;
}

void
dict_report(dict_t * d)
{
    E_INFO_NOFN("Initialization of dict_t, report:\n");
    E_INFO_NOFN("Max word: %d\n", d->max_words);
    E_INFO_NOFN("No of word: %d\n", d->n_word);
    E_INFO_NOFN("\n");
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

#ifndef _S3_DICT_H_
#define _S3_DICT_H_

/** \file dict.h
 * \brief Operations on dictionary. 
 */

/* SphinxBase headers. */
#include <sphinxbase/hash_table.h>

/* Local headers. */
#include "s3types.h"
#include "bin_mdef.h"
#include "pocketsphinx_export.h"

#define S3DICT_INC_SZ 4096

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

/** 
    \struct dictword_t
    \brief a structure for one dictionary word. 
*/
typedef struct {
    char *word;		/**< Ascii word string */
    s3cipid_t *ciphone;	/**< Pronunciation */
    int32 pronlen;	/**< Pronunciation length */
    s3wid_t alt;	/**< Next alternative pronunciation id, NOT_S3WID if none */
    s3wid_t basewid;	/**< Base pronunciation id */
} dictword_t;

/** 
    \struct dict_t
    \brief a structure for a dictionary. 
*/

typedef struct {
    int refcnt;
    bin_mdef_t *mdef;	/**< Model definition used for phone IDs; NULL if none used */
    dictword_t *word;	/**< Array of entries in dictionary */
    hash_table_t *ht;	/**< Hash table for mapping word strings to word ids */
    int32 max_words;	/**< #Entries allocated in dict, including empty slots */
    int32 n_word;	/**< #Occupied entries in dict; ie, excluding empty slots */
    int32 filler_start;	/**< First filler word id (read from filler dict) */
    int32 filler_end;	/**< Last filler word id (read from filler dict) */
    s3wid_t startwid;	/**< FOR INTERNAL-USE ONLY */
    s3wid_t finishwid;	/**< FOR INTERNAL-USE ONLY */
    s3wid_t silwid;	/**< FOR INTERNAL-USE ONLY */
    int nocase;
} dict_t;


/**
 * Initialize a new dictionary.
 *
 * If config and mdef are supplied, then the dictionary will be read
 * from the files specified by the -dict and -fdict options in config,
 * with case sensitivity determined by the -dictcase option.
 *
 * Otherwise an empty case-sensitive dictionary will be created.
 *
 * Return ptr to dict_t if successful, NULL otherwise.
 */
dict_t *dict_init(cmd_ln_t *config, /**< Configuration (-dict, -fdict, -dictcase) or NULL */
                  bin_mdef_t *mdef  /**< For looking up CI phone IDs (or NULL) */
    );

/**
 * Write dictionary to a file.
 */
int dict_write(dict_t *dict, char const *filename, char const *format);

/** Return word id for given word string if present.  Otherwise return BAD_S3WID */
POCKETSPHINX_EXPORT
s3wid_t dict_wordid(dict_t *d, const char *word);

/**
 * Return 1 if w is a filler word, 0 if not.  A filler word is one that was read in from the
 * filler dictionary; however, sentence START and FINISH words are not filler words.
 */
int dict_filler_word(dict_t *d,  /**< The dictionary structure */
                     s3wid_t w     /**< The word ID */
    );

/**
 * Test if w is a "real" word, i.e. neither a filler word nor START/FINISH.
 */
POCKETSPHINX_EXPORT
int dict_real_word(dict_t *d,  /**< The dictionary structure */
                   s3wid_t w     /**< The word ID */
    );

/**
 * Add a word with the given ciphone pronunciation list to the dictionary.
 * Return value: Result word id if successful, BAD_S3WID otherwise
 */
s3wid_t dict_add_word(dict_t *d,          /**< The dictionary structure. */
                      char const *word,   /**< The word. */
                      s3cipid_t const *p, /**< The pronunciation. */
                      int32 np            /**< Number of phones. */
    );

/**
 * Return value: CI phone string for the given word, phone position.
 */
const char *dict_ciphone_str(dict_t *d,	/**< In: Dictionary to look up */
                             s3wid_t wid,	/**< In: Component word being looked up */
                             int32 pos   	/**< In: Pronunciation phone position */
    );

/** Packaged macro access to dictionary members */
#define dict_size(d)		((d)->n_word)
#define dict_num_fillers(d)   (dict_filler_end(d) - dict_filler_start(d))
/**
 * Number of "real words" in the dictionary.
 *
 * This is the number of words that are not fillers, <s>, or </s>.
 */
#define dict_num_real_words(d)                                          \
    (dict_size(d) - (dict_filler_end(d) - dict_filler_start(d)) - 2)
#define dict_basewid(d,w)	((d)->word[w].basewid)
#define dict_wordstr(d,w)	((w) < 0 ? NULL : (d)->word[w].word)
#define dict_basestr(d,w)	((d)->word[dict_basewid(d,w)].word)
#define dict_nextalt(d,w)	((d)->word[w].alt)
#define dict_pronlen(d,w)	((d)->word[w].pronlen) 
#define dict_pron(d,w,p)	((d)->word[w].ciphone[p]) /**< The CI phones of the word w at position p */
#define dict_filler_start(d)	((d)->filler_start)
#define dict_filler_end(d)	((d)->filler_end)
#define dict_startwid(d)	((d)->startwid)
#define dict_finishwid(d)	((d)->finishwid)
#define dict_silwid(d)		((d)->silwid)
#define dict_is_single_phone(d,w)	((d)->word[w].pronlen == 1)
#define dict_first_phone(d,w)	((d)->word[w].ciphone[0])
#define dict_second_phone(d,w)	((d)->word[w].ciphone[1])
#define dict_second_last_phone(d,w)	((d)->word[w].ciphone[(d)->word[w].pronlen - 2])
#define dict_last_phone(d,w)	((d)->word[w].ciphone[(d)->word[w].pronlen - 1])

/* Hard-coded special words */
#define S3_START_WORD		"<s>"
#define S3_FINISH_WORD		"</s>"
#define S3_SILENCE_WORD		"<sil>"
#define S3_UNKNOWN_WORD		"<UNK>"

/**
 * If the given word contains a trailing "(....)" (i.e., a Sphinx-II style alternative
 * pronunciation specification), strip that trailing portion from it.  Note that the given
 * string is modified.
 * Return value: If string was modified, the character position at which the original string
 * was truncated; otherwise -1.
 */
int32 dict_word2basestr(char *word);

/**
 * Retain a pointer to an dict_t.
 */
dict_t *dict_retain(dict_t *d);

/**
 * Release a pointer to a dictionary.
 */
int dict_free(dict_t *d);

/** Report a dictionary structure */
void dict_report(dict_t *d /**< A dictionary structure */
    );

#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif

#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * fsg_history.c -- FSG Viterbi decode history
 * 
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1999 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * 
 * 25-Feb-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Started..
 */

/* System headers. */
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>

/* Local headers. */
#include "fsg_search_internal.h"
#include "fsg_history.h"


#define __FSG_DBG__	0


fsg_history_t *
fsg_history_init(fsg_model_t * fsg, dict_t *dict)
{
    fsg_history_t *h;

    h = (fsg_history_t *) ckd_calloc(1, sizeof(fsg_history_t));
    h->fsg = fsg;
    h->entries = blkarray_list_init();

    if (fsg && dict) {
        h->n_ciphone = bin_mdef_n_ciphone(dict->mdef);
        h->frame_entries =
            (glist_t **) ckd_calloc_2d(fsg_model_n_state(fsg),
                                       bin_mdef_n_ciphone(dict->mdef),
                                       sizeof(**h->frame_entries));
    }
    else {
        h->frame_entries = NULL;
    }

    return h;
}

void
fsg_history_free(fsg_history_t *h)
{
    int32 s, lc, ns, np;
    gnode_t *gn;

    if (h->fsg) {
        ns = fsg_model_n_state(h->fsg);
        np = h->n_ciphone;

        for (s = 0; s < ns; s++) {
            for (lc = 0; lc < np; lc++) {
                for (gn = h->frame_entries[s][lc]; gn; gn = gnode_next(gn)) {
                    ckd_free(gnode_ptr(gn));
                }
                glist_free(h->frame_entries[s][lc]);
            }
        }
    }
    ckd_free_2d(h->frame_entries);
    blkarray_list_free(h->entries);
    ckd_free(h);
}


void
fsg_history_set_fsg(fsg_history_t *h, fsg_model_t *fsg, dict_t *dict)
{
    if (blkarray_list_n_valid(h->entries) != 0) {
        E_WARN("Switching FSG while history not empty; history cleared\n");
        blkarray_list_reset(h->entries);
    }

    if (h->frame_entries)
        ckd_free_2d((void **) h->frame_entries);
    h->frame_entries = NULL;
    h->fsg = fsg;

    if (fsg && dict) {
        h->n_ciphone = bin_mdef_n_ciphone(dict->mdef);
        h->frame_entries =
            (glist_t **) ckd_calloc_2d(fsg_model_n_state(fsg),
                                       bin_mdef_n_ciphone(dict->mdef),
                                       sizeof(glist_t));
    }
}


void
fsg_history_entry_add(fsg_history_t * h,
                      fsg_link_t * link,
                      int32 frame, int32 score, int32 pred,
                      int32 lc, fsg_pnode_ctxt_t rc)
{
    fsg_hist_entry_t *entry, *new_entry;
    int32 s;
    gnode_t *gn, *prev_gn;

    /* Skip the optimization for the initial dummy entries; always enter them */
    if (frame < 0) {
        new_entry =
            (fsg_hist_entry_t *) ckd_calloc(1, sizeof(fsg_hist_entry_t));
        new_entry->fsglink = link;
        new_entry->frame = frame;
        new_entry->score = score;
        new_entry->pred = pred;
        new_entry->lc = lc;
        new_entry->rc = rc;

        blkarray_list_append(h->entries, (void *) new_entry);
        return;
    }

    s = fsg_link_to_state(link);

    /* Locate where this entry should be inserted in frame_entries[s][lc] */
    prev_gn = NULL;
    for (gn = h->frame_entries[s][lc]; gn; gn = gnode_next(gn)) {
        entry = (fsg_hist_entry_t *) gnode_ptr(gn);

        if (score BETTER_THAN entry->score)
            break;              /* Found where to insert new entry */

        /* Existing entry score not worse than new score */
        if (FSG_PNODE_CTXT_SUB(&rc, &(entry->rc)) == 0)
            return;             /* rc set reduced to 0; new entry can be ignored */

        prev_gn = gn;
    }

    /* Create new entry after prev_gn (if prev_gn is NULL, at head) */
    new_entry =
        (fsg_hist_entry_t *) ckd_calloc(1, sizeof(fsg_hist_entry_t));
    new_entry->fsglink = link;
    new_entry->frame = frame;
    new_entry->score = score;
    new_entry->pred = pred;
    new_entry->lc = lc;
    new_entry->rc = rc;         /* Note: rc set must be non-empty at this point */

    if (!prev_gn) {
        h->frame_entries[s][lc] = glist_add_ptr(h->frame_entries[s][lc],
                                                (void *) new_entry);
        prev_gn = h->frame_entries[s][lc];
    }
    else
        prev_gn = glist_insert_ptr(prev_gn, (void *) new_entry);

    /*
     * Update the rc set of all the remaining entries in the list.  At this
     * point, gn is the entry, if any, immediately following new entry.
     */
    while (gn) {
        entry = (fsg_hist_entry_t *) gnode_ptr(gn);

        if (FSG_PNODE_CTXT_SUB(&(entry->rc), &rc) == 0) {
            /* rc set of entry reduced to 0; can prune this entry */
            ckd_free((void *) entry);
            gn = gnode_free(gn, prev_gn);
        }
        else {
            prev_gn = gn;
            gn = gnode_next(gn);
        }
    }
}


/*
 * Transfer the surviving history entries for this frame into the permanent
 * history table.
 */
void
fsg_history_end_frame(fsg_history_t * h)
{
    int32 s, lc, ns, np;
    gnode_t *gn;
    fsg_hist_entry_t *entry;

    ns = fsg_model_n_state(h->fsg);
    np = h->n_ciphone;

    for (s = 0; s < ns; s++) {
        for (lc = 0; lc < np; lc++) {
            for (gn = h->frame_entries[s][lc]; gn; gn = gnode_next(gn)) {
                entry = (fsg_hist_entry_t *) gnode_ptr(gn);
                blkarray_list_append(h->entries, (void *) entry);
            }

            glist_free(h->frame_entries[s][lc]);
            h->frame_entries[s][lc] = NULL;
        }
    }
}


fsg_hist_entry_t *
fsg_history_entry_get(fsg_history_t * h, int32 id)
{
    blkarray_list_t *entries;
    int32 r, c;

    entries = h->entries;

    if (id >= blkarray_list_n_valid(entries))
        return NULL;

    r = id / blkarray_list_blksize(entries);
    c = id - (r * blkarray_list_blksize(entries));

    return ((fsg_hist_entry_t *) blkarray_list_ptr(entries, r, c));
}


void
fsg_history_reset(fsg_history_t * h)
{
    blkarray_list_reset(h->entries);
}


int32
fsg_history_n_entries(fsg_history_t * h)
{
    return (blkarray_list_n_valid(h->entries));
}

void
fsg_history_utt_start(fsg_history_t * h)
{
    int32 s, lc, ns, np;

    assert(blkarray_list_n_valid(h->entries) == 0);
    assert(h->frame_entries);

    ns = fsg_model_n_state(h->fsg);
    np = h->n_ciphone;

    for (s = 0; s < ns; s++) {
        for (lc = 0; lc < np; lc++) {
            assert(h->frame_entries[s][lc] == NULL);
        }
    }
}

void
fsg_history_utt_end(fsg_history_t * h)
{
}

void
fsg_history_print(fsg_history_t *h, dict_t *dict) 
{
    int bpidx, bp;
    
    for (bpidx = 0; bpidx < blkarray_list_n_valid(h->entries); bpidx++) {
        bp = bpidx;
        printf("History entry: ");
        while (bp > 0) {
            fsg_hist_entry_t *hist_entry = fsg_history_entry_get(h, bp);
	    fsg_link_t *fl = fsg_hist_entry_fsglink(hist_entry);
    	    char const *baseword;
    	    int32 wid;
    	    bp = fsg_hist_entry_pred(hist_entry);
    	    wid = fsg_link_wid(fl);

    	    if (fl == NULL)
        	    continue;

    	    baseword = fsg_model_word_str(h->fsg, wid);

    	    printf("%s(%d->%d:%d) ", baseword, 
    				     fsg_link_from_state(hist_entry->fsglink), 
    				     fsg_link_to_state(hist_entry->fsglink), 
    				     hist_entry->frame);
	}
	printf("\n");
    }
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * fsg_history.h -- FSG Viterbi decode history
 * 
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1999 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * 
 * $Log: fsg_history.h,v $
 * Revision 1.1.1.1  2006/05/23 18:45:02  dhuggins
 * re-importation
 *
 * Revision 1.1  2004/07/16 00:57:12  egouvea
 * Added Ravi's implementation of FSG support.
 *
 * Revision 1.7  2004/07/07 22:30:35  rkm
 * *** empty log message ***
 *
 * Revision 1.6  2004/07/07 13:56:33  rkm
 * Added reporting of (acoustic score - best senone score)/frame
 *
 * Revision 1.5  2004/06/25 14:49:08  rkm
 * Optimized size of history table and speed of word transitions by maintaining only best scoring word exits at each state
 *
 * Revision 1.4  2004/06/23 20:32:16  rkm
 * *** empty log message ***
 *
 * Revision 1.3  2004/05/27 15:16:08  rkm
 * *** empty log message ***
 *
 * 
 * 25-Feb-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Started, based on S3.3 version.
 */


#ifndef __S2_FSG_HISTORY_H__
#define __S2_FSG_HISTORY_H__


/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>
#include <sphinxbase/fsg_model.h>

/* Local headers. */
#include "blkarray_list.h"
#include "fsg_lextree.h"
#include "dict.h"

/*
 * The Viterbi history structure.  This is a tree, with the root at the
 * FSG start state, at frame 0, with a null predecessor.
 */

/*
 * A single Viterbi history entry
 */
typedef struct fsg_hist_entry_s {
    fsg_link_t *fsglink;		/* Link taken result in this entry */
    int32 score;			/* Total path score at the end of this
                                           transition */
    int32 pred; 			/* Predecessor entry; -1 if none */
    frame_idx_t frame;			/* Ending frame for this entry */
    int16 lc;			        /* Left context provided by this entry to
				           succeeding words */
    fsg_pnode_ctxt_t rc;		/* Possible right contexts to which this entry
                                           applies */
} fsg_hist_entry_t;

/* Access macros */
#define fsg_hist_entry_fsglink(v)	((v)->fsglink)
#define fsg_hist_entry_frame(v)		((v)->frame)
#define fsg_hist_entry_score(v)		((v)->score)
#define fsg_hist_entry_pred(v)		((v)->pred)
#define fsg_hist_entry_lc(v)		((v)->lc)
#define fsg_hist_entry_rc(v)		((v)->rc)


/*
 * The entire tree of history entries (fsg_history_t.entries).
 * Optimization: In a given frame, there may be several history entries, with
 * the same left and right phonetic context, terminating in a particular state.
 * Only the best scoring one of these needs to be saved, since everything else
 * will be pruned according to the Viterbi algorithm.  frame_entries is used
 * temporarily in each frame to determine these best scoring entries in that
 * frame.  Only the ones not pruned are transferred to entries at the end of
 * the frame.  However, null transitions are a problem since they create
 * entries that depend on entries created in the CURRENT frame.  Hence, this
 * pruning is done in two stages: first for the non-null transitions, and then
 * for the null transitions alone.  (This solution is sub-optimal, and can be
 * improved with a little more work.  SMOP.)
 * Why is frame_entries a list?  Each entry has a unique terminating state,
 * and has a unique lc CIphone.  But it has a SET of rc CIphones.
 * frame_entries[s][lc] is an ordered list of entries created in the current
 * frame, terminating in state s, and with left context lc.  The list is in
 * descending order of path score.  When a new entry with (s,lc) arrives,
 * its position in the list is determined.  Then its rc set is modified by
 * subtracting the union of the rc's of all its predecessors (i.e., better
 * scoring entries).  If the resulting rc set is empty, the entry is discarded.
 * Otherwise, it is inserted, and the rc sets of all downstream entries in the
 * list are updated by subtracting the new entry's rc.  If any of them becomes
 * empty, it is also discarded.
 * As mentioned earlier, this procedure is applied in two stages, for the
 * non-null transitions, and the null transitions, separately.
 */
typedef struct fsg_history_s {
    fsg_model_t *fsg;		/* The FSG for which this object applies */
    blkarray_list_t *entries;	/* A list of history table entries; the root
				   entry is the first element of the list */
    glist_t **frame_entries;
    int n_ciphone;
} fsg_history_t;


/*
 * One-time intialization: Allocate and return an initially empty history
 * module.
 */
fsg_history_t *fsg_history_init(fsg_model_t *fsg, dict_t *dict);

void fsg_history_utt_start(fsg_history_t *h);

void fsg_history_utt_end(fsg_history_t *h);


/*
 * Create a history entry recording the completion of the given FSG
 * transition, at the end of the given frame, with the given score, and
 * the given predecessor history entry.
 * The entry is initially temporary, and may be superseded by another
 * with a higher score.  The surviving entries must be transferred to
 * the main history table, via fsg_history_end_frame().
 */
void fsg_history_entry_add (fsg_history_t *h,
			    fsg_link_t *l,	/* FSG transition */
			    int32 frame,
			    int32 score,
			    int32 pred,
			    int32 lc,
			    fsg_pnode_ctxt_t rc);

/*
 * Transfer the surviving history entries for this frame into the permanent
 * history table.  This function can be called several times during a frame.
 * Each time, the entries surviving so far are transferred, and the temporary
 * lists cleared.  This feature is used to handle the entries due to non-null
 * transitions and null transitions separately.
 */
void fsg_history_end_frame (fsg_history_t *h);


/* Clear the hitory table */
void fsg_history_reset (fsg_history_t *h);


/* Return the number of valid entries in the given history table */
int32 fsg_history_n_entries (fsg_history_t *h);

/*
 * Return a ptr to the history entry for the given ID; NULL if there is no
 * such entry.
 */
fsg_hist_entry_t *fsg_history_entry_get(fsg_history_t *h, int32 id);


/*
 * Switch the FSG associated with the given history module.  Should be done
 * when the history list is empty.  If not empty, the list is cleared.
 */
void fsg_history_set_fsg (fsg_history_t *h, fsg_model_t *fsg, dict_t *dict);

/* Free the given Viterbi search history object */
void fsg_history_free (fsg_history_t *h);

/* Print the entire history */
void fsg_history_print(fsg_history_t *h, dict_t *dict);
				     
#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/**
 * @file fsg_lextree.c
 * @brief The collection of all the lextrees for the entire FSM.
 * @author M K Ravishankar <rkm@cs.cmu.edu>
 * @author Bhiksha Raj <bhiksha@cs.cmu.edu>
 */

/* System headers. */
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "fsg_lextree.h"

#define __FSG_DBG__		0

/* A linklist structure that is actually used to build local lextrees at grammar nodes */
typedef struct fsg_glist_linklist_t {
    int32    ci, rc;
    glist_t  glist;
    struct   fsg_glist_linklist_t *next;
} fsg_glist_linklist_t;

/**
 * Build the phone lextree for all transitions out of state from_state.
 * Return the root node of this tree.
 * Also, return a linear linked list of all allocated fsg_pnode_t nodes in
 * *alloc_head (for memory management purposes).
 */
static fsg_pnode_t *fsg_psubtree_init(fsg_lextree_t *tree,
                                      fsg_model_t *fsg,
                                      int32 from_state,
                                      fsg_pnode_t **alloc_head);

/**
 * Free the given lextree.  alloc_head: head of linear list of allocated
 * nodes updated by fsg_psubtree_init().
 */
static void fsg_psubtree_free(fsg_pnode_t *alloc_head);

/**
 * Dump the list of nodes in the given lextree to the given file.  alloc_head:
 * head of linear list of allocated nodes updated by fsg_psubtree_init().
 */
static void fsg_psubtree_dump(fsg_lextree_t *tree, fsg_pnode_t *root, FILE *fp);

/**
 * Compute the left and right context CIphone sets for each state.
 */
static void
fsg_lextree_lc_rc(fsg_lextree_t *lextree)
{
    int32 s, i, j;
    int32 n_ci;
    fsg_model_t *fsg;
    int32 silcipid;
    int32 len;

    silcipid = bin_mdef_silphone(lextree->mdef);
    assert(silcipid >= 0);
    n_ci = bin_mdef_n_ciphone(lextree->mdef);

    fsg = lextree->fsg;
    /*
     * lextree->lc[s] = set of left context CIphones for state s.  Similarly, rc[s]
     * for right context CIphones.
     */
    lextree->lc = ckd_calloc_2d(fsg->n_state, n_ci + 1, sizeof(**lextree->lc));
    lextree->rc = ckd_calloc_2d(fsg->n_state, n_ci + 1, sizeof(**lextree->rc));
    E_INFO("Allocated %d bytes (%d KiB) for left and right context phones\n",
           fsg->n_state * (n_ci + 1) * 2,
           fsg->n_state * (n_ci + 1) * 2 / 1024);


    for (s = 0; s < fsg->n_state; s++) {
        fsg_arciter_t *itor;
        for (itor = fsg_model_arcs(fsg, s); itor; itor = fsg_arciter_next(itor)) {
            fsg_link_t *l = fsg_arciter_get(itor);
            int32 dictwid; /**< Dictionary (not FSG) word ID!! */

            if (fsg_link_wid(l) >= 0) {
                dictwid = dict_wordid(lextree->dict,
                                      fsg_model_word_str(lextree->fsg, l->wid));

                /*
                 * Add the first CIphone of l->wid to the rclist of state s, and
                 * the last CIphone to lclist of state d.
                 * (Filler phones are a pain to deal with.  There is no direct
                 * marking of a filler phone; but only filler words are supposed to
                 * use such phones, so we use that fact.  HACK!!  FRAGILE!!)
                 */
                if (fsg_model_is_filler(fsg, fsg_link_wid(l))) {
                    /* Filler phone; use silence phone as context */
                    lextree->rc[fsg_link_from_state(l)][silcipid] = 1;
                    lextree->lc[fsg_link_to_state(l)][silcipid] = 1;
                }
                else {
                    len = dict_pronlen(lextree->dict, dictwid);
                    lextree->rc[fsg_link_from_state(l)][dict_pron(lextree->dict, dictwid, 0)] = 1;
                    lextree->lc[fsg_link_to_state(l)][dict_pron(lextree->dict, dictwid, len - 1)] = 1;
                }
            }
        }
    }

    for (s = 0; s < fsg->n_state; s++) {
        /*
         * Add SIL phone to the lclist and rclist of each state.  Strictly
         * speaking, only needed at start and final states, respectively, but
         * all states considered since the user may change the start and final
         * states.  In any case, most applications would have a silence self
         * loop at each state, hence these would be needed anyway.
         */
        lextree->lc[s][silcipid] = 1;
        lextree->rc[s][silcipid] = 1;
    }

    /*
     * Propagate lc and rc lists past null transitions.  (Since FSG contains
     * null transitions closure, no need to worry about a chain of successive
     * null transitions.  Right??)
     *
     * This can't be joined with the previous loop because we first calculate 
     * contexts and only then we can propagate them.
     */
    for (s = 0; s < fsg->n_state; s++) {
        fsg_arciter_t *itor;
        for (itor = fsg_model_arcs(fsg, s); itor; itor = fsg_arciter_next(itor)) {
            fsg_link_t *l = fsg_arciter_get(itor);
            if (fsg_link_wid(l) < 0) {

                /*
                 * lclist(d) |= lclist(s), because all the words ending up at s, can
                 * now also end at d, becoming the left context for words leaving d.
                 */
                for (i = 0; i < n_ci; i++)
                    lextree->lc[fsg_link_to_state(l)][i] |= lextree->lc[fsg_link_from_state(l)][i];
                /*
                 * Similarly, rclist(s) |= rclist(d), because all the words leaving d
                 * can equivalently leave s, becoming the right context for words
                 * ending up at s.
                 */
                for (i = 0; i < n_ci; i++)
                    lextree->rc[fsg_link_from_state(l)][i] |= lextree->rc[fsg_link_to_state(l)][i];
            }
        }
    }

    /* Convert the bit-vector representation into a list */
    for (s = 0; s < fsg->n_state; s++) {
        j = 0;
        for (i = 0; i < n_ci; i++) {
            if (lextree->lc[s][i]) {
                lextree->lc[s][j] = i;
                j++;
            }
        }
        lextree->lc[s][j] = -1;     /* Terminate the list */

        j = 0;
        for (i = 0; i < n_ci; i++) {
            if (lextree->rc[s][i]) {
                lextree->rc[s][j] = i;
                j++;
            }
        }
        lextree->rc[s][j] = -1;     /* Terminate the list */
    }
}

/*
 * For now, allocate the entire lextree statically.
 */
fsg_lextree_t *
fsg_lextree_init(fsg_model_t * fsg, dict_t *dict, dict2pid_t *d2p,
                 bin_mdef_t *mdef, hmm_context_t *ctx,
                 int32 wip, int32 pip)
{
    int32 s, n_leaves;
    fsg_lextree_t *lextree;
    fsg_pnode_t *pn;

    lextree = ckd_calloc(1, sizeof(fsg_lextree_t));
    lextree->fsg = fsg;
    lextree->root = ckd_calloc(fsg_model_n_state(fsg),
                               sizeof(fsg_pnode_t *));
    lextree->alloc_head = ckd_calloc(fsg_model_n_state(fsg),
                                     sizeof(fsg_pnode_t *));
    lextree->ctx = ctx;
    lextree->dict = dict;
    lextree->d2p = d2p;
    lextree->mdef = mdef;
    lextree->wip = wip;
    lextree->pip = pip;

    /* Compute lc and rc for fsg. */
    fsg_lextree_lc_rc(lextree);

    /* Create lextree for each state, i.e. an HMM network that
     * represents words for all arcs exiting that state.  Note that
     * for a dense grammar such as an N-gram model, this will
     * rapidly exhaust all available memory. */
    lextree->n_pnode = 0;
    n_leaves = 0;
    for (s = 0; s < fsg_model_n_state(fsg); s++) {
        lextree->root[s] =
            fsg_psubtree_init(lextree, fsg, s, &(lextree->alloc_head[s]));

        for (pn = lextree->alloc_head[s]; pn; pn = pn->alloc_next) {
            lextree->n_pnode++;
            if (pn->leaf)
                ++n_leaves;
        }
    }
    E_INFO("%d HMM nodes in lextree (%d leaves)\n",
           lextree->n_pnode, n_leaves);
    E_INFO("Allocated %d bytes (%d KiB) for all lextree nodes\n",
           lextree->n_pnode * sizeof(fsg_pnode_t),
           lextree->n_pnode * sizeof(fsg_pnode_t) / 1024);
    E_INFO("Allocated %d bytes (%d KiB) for lextree leafnodes\n",
           n_leaves * sizeof(fsg_pnode_t),
           n_leaves * sizeof(fsg_pnode_t) / 1024);

#if __FSG_DBG__
    fsg_lextree_dump(lextree, stdout);
#endif

    return lextree;
}


void
fsg_lextree_dump(fsg_lextree_t * lextree, FILE * fp)
{
    int32 s;

    for (s = 0; s < fsg_model_n_state(lextree->fsg); s++) {
        fprintf(fp, "State %5d root %p\n", s, lextree->root[s]);
        fsg_psubtree_dump(lextree, lextree->root[s], fp);
    }
    fflush(fp);
}


void
fsg_lextree_free(fsg_lextree_t * lextree)
{
    int32 s;

    if (lextree == NULL)
        return;

    if (lextree->fsg)
        for (s = 0; s < fsg_model_n_state(lextree->fsg); s++)
            fsg_psubtree_free(lextree->alloc_head[s]);

    ckd_free_2d(lextree->lc);
    ckd_free_2d(lextree->rc);
    ckd_free(lextree->root);
    ckd_free(lextree->alloc_head);
    ckd_free(lextree);
}

/******************************
 * psubtree stuff starts here *
 ******************************/

void fsg_glist_linklist_free(fsg_glist_linklist_t *glist)
{
    if (glist) {
        fsg_glist_linklist_t *nxtglist;
        if (glist->glist)
            glist_free(glist->glist);
        nxtglist = glist->next;
        while (nxtglist) {
            ckd_free(glist);
            glist = nxtglist;
            if (glist->glist)
                glist_free(glist->glist);
            nxtglist = glist->next;
        }
        ckd_free(glist);
    }
    return;
}

void
fsg_pnode_add_all_ctxt(fsg_pnode_ctxt_t * ctxt)
{
    int32 i;

    for (i = 0; i < FSG_PNODE_CTXT_BVSZ; i++)
        ctxt->bv[i] = 0xffffffff;
}

uint32 fsg_pnode_ctxt_sub_generic(fsg_pnode_ctxt_t *src, fsg_pnode_ctxt_t *sub)
{
    int32 i;
    uint32 res = 0;
    
    for (i = 0; i < FSG_PNODE_CTXT_BVSZ; i++)
        res |= (src->bv[i] = ~(sub->bv[i]) & src->bv[i]);
    return res;
}


/*
 * fsg_pnode_ctxt_sub(fsg_pnode_ctxt_t * src, fsg_pnode_ctxt_t * sub)
 * This has been moved into a macro in fsg_psubtree.h 
 * because it is called so frequently!
 */


/*
 * Add the word emitted by the given transition (fsglink) to the given lextree
 * (rooted at root), and return the new lextree root.  (There may actually be
 * several root nodes, maintained in a linked list via fsg_pnode_t.sibling.
 * "root" is the head of this list.)
 * lclist, rclist: sets of left and right context phones for this link.
 * alloc_head: head of a linear list of all allocated pnodes for the parent
 * FSG state, kept elsewhere and updated by this routine.
 */
static fsg_pnode_t *
psubtree_add_trans(fsg_lextree_t *lextree, 
                   fsg_pnode_t * root,
                   fsg_glist_linklist_t **curglist,
                   fsg_link_t * fsglink,
                   int16 *lclist, int16 *rclist,
                   fsg_pnode_t ** alloc_head)
{
    int32 silcipid;             /* Silence CI phone ID */
    int32 pronlen;              /* Pronunciation length */
    int32 wid;                  /* FSG (not dictionary!!) word ID */
    int32 dictwid;              /* Dictionary (not FSG!!) word ID */
    int32 ssid;                 /* Senone Sequence ID */
    int32 tmatid;
    gnode_t *gn;
    fsg_pnode_t *pnode, *pred, *head;
    int32 n_ci, p, lc, rc;
    glist_t lc_pnodelist;       /* Temp pnodes list for different left contexts */
    glist_t rc_pnodelist;       /* Temp pnodes list for different right contexts */
    int32 i, j;
    int n_lc_alloc = 0, n_int_alloc = 0, n_rc_alloc = 0;

    silcipid = bin_mdef_silphone(lextree->mdef);
    n_ci = bin_mdef_n_ciphone(lextree->mdef);

    wid = fsg_link_wid(fsglink);
    assert(wid >= 0);           /* Cannot be a null transition */
    dictwid = dict_wordid(lextree->dict,
                          fsg_model_word_str(lextree->fsg, wid));
    pronlen = dict_pronlen(lextree->dict, dictwid);
    assert(pronlen >= 1);

    assert(lclist[0] >= 0);     /* At least one phonetic context provided */
    assert(rclist[0] >= 0);

    head = *alloc_head;
    pred = NULL;

    if (pronlen == 1) {         /* Single-phone word */
        int ci = dict_first_phone(lextree->dict, dictwid);
        /* Only non-filler words are mpx */
        if (dict_filler_word(lextree->dict, dictwid)) {
            /*
             * Left diphone ID for single-phone words already assumes SIL is right
             * context; only left contexts need to be handled.
             */
            lc_pnodelist = NULL;

            for (i = 0; lclist[i] >= 0; i++) {
                lc = lclist[i];
                ssid = dict2pid_lrdiph_rc(lextree->d2p, ci, lc, silcipid);
                tmatid = bin_mdef_pid2tmatid(lextree->mdef, dict_first_phone(lextree->dict, dictwid));
                /* Check if this ssid already allocated for some other context */
                for (gn = lc_pnodelist; gn; gn = gnode_next(gn)) {
                    pnode = (fsg_pnode_t *) gnode_ptr(gn);

                    if (hmm_nonmpx_ssid(&pnode->hmm) == ssid) {
                        /* already allocated; share it for this context phone */
                        fsg_pnode_add_ctxt(pnode, lc);
                        break;
                    }
                }

                if (!gn) {      /* ssid not already allocated */
                    pnode =
                        (fsg_pnode_t *) ckd_calloc(1, sizeof(fsg_pnode_t));
                    pnode->ctx = lextree->ctx;
                    pnode->next.fsglink = fsglink;
                    pnode->logs2prob =
                        (fsg_link_logs2prob(fsglink) >> SENSCR_SHIFT)
                        + lextree->wip + lextree->pip;
                    pnode->ci_ext = dict_first_phone(lextree->dict, dictwid);
                    pnode->ppos = 0;
                    pnode->leaf = TRUE;
                    pnode->sibling = root;      /* All root nodes linked together */
                    fsg_pnode_add_ctxt(pnode, lc);      /* Initially zeroed by calloc above */
                    pnode->alloc_next = head;
                    head = pnode;
                    root = pnode;
                    ++n_lc_alloc;

                    hmm_init(lextree->ctx, &pnode->hmm, FALSE, ssid, tmatid);

                    lc_pnodelist =
                        glist_add_ptr(lc_pnodelist, (void *) pnode);
                }
            }

            glist_free(lc_pnodelist);
        }
        else {                  /* Filler word; no context modelled */
            ssid = bin_mdef_pid2ssid(lextree->mdef, ci); /* probably the same... */
            tmatid = bin_mdef_pid2tmatid(lextree->mdef, ci);

            pnode = (fsg_pnode_t *) ckd_calloc(1, sizeof(fsg_pnode_t));
            pnode->ctx = lextree->ctx;
            pnode->next.fsglink = fsglink;
            pnode->logs2prob = (fsg_link_logs2prob(fsglink) >> SENSCR_SHIFT)
                + lextree->wip + lextree->pip;
            pnode->ci_ext = silcipid;   /* Presents SIL as context to neighbors */
            pnode->ppos = 0;
            pnode->leaf = TRUE;
            pnode->sibling = root;
            fsg_pnode_add_all_ctxt(&(pnode->ctxt));
            pnode->alloc_next = head;
            head = pnode;
            root = pnode;
            ++n_int_alloc;

            hmm_init(lextree->ctx, &pnode->hmm, FALSE, ssid, tmatid);
        }
    }
    else {                      /* Multi-phone word */
        fsg_pnode_t **ssid_pnode_map;       /* Temp array of ssid->pnode mapping */
        ssid_pnode_map =
            (fsg_pnode_t **) ckd_calloc(n_ci, sizeof(fsg_pnode_t *));
        lc_pnodelist = NULL;
        rc_pnodelist = NULL;

        for (p = 0; p < pronlen; p++) {
            int ci = dict_pron(lextree->dict, dictwid, p);
            if (p == 0) {       /* Root phone, handle required left contexts */
                /* Find if we already have an lc_pnodelist for the first phone of this word */
		fsg_glist_linklist_t *glist;

                rc = dict_pron(lextree->dict, dictwid, 1);
		for (glist = *curglist;
                     glist && glist->glist;
                     glist = glist->next) {
		    if (glist->ci == ci && glist->rc == rc)
			break;
		}
		if (glist && glist->glist) {
		    assert(glist->ci == ci && glist->rc == rc);
		    /* We've found a valid glist. Hook to it and move to next phoneme */
                    E_DEBUG(2,("Found match for (%d,%d)\n", ci, rc));
		    lc_pnodelist = glist->glist;
                    /* Set the predecessor node for the future tree first */
		    pred = (fsg_pnode_t *) gnode_ptr(lc_pnodelist);
		    continue;
		}
		else {
		    /* Two cases that can bring us here
		     * a. glist == NULL, i.e. end of current list. Create new entry.
		     * b. glist->glist == NULL, i.e. first entry into list.
		     */
		    if (glist == NULL) { /* Case a; reduce it to case b by allocing glist */
		        glist = (fsg_glist_linklist_t*) ckd_calloc(1, sizeof(fsg_glist_linklist_t));
			glist->next = *curglist;
                        *curglist = glist;
		    }
		    glist->ci = ci;
                    glist->rc = rc;
		    lc_pnodelist = glist->glist = NULL; /* Gets created below */
		}

                for (i = 0; lclist[i] >= 0; i++) {
                    lc = lclist[i];
                    ssid = dict2pid_ldiph_lc(lextree->d2p, ci, rc, lc);
                    tmatid = bin_mdef_pid2tmatid(lextree->mdef, dict_first_phone(lextree->dict, dictwid));
                    /* Compression is not done by d2p, so we do it
                     * here.  This might be slow, but it might not
                     * be... we'll see. */
                    pnode = ssid_pnode_map[0];
                    for (j = 0; j < n_ci && ssid_pnode_map[j] != NULL; ++j) {
                        pnode = ssid_pnode_map[j];
                        if (hmm_nonmpx_ssid(&pnode->hmm) == ssid)
                            break;
                    }
                    assert(j < n_ci);
                    if (!pnode) {       /* Allocate pnode for this new ssid */
                        pnode =
                            (fsg_pnode_t *) ckd_calloc(1,
                                                       sizeof
                                                       (fsg_pnode_t));
                        pnode->ctx = lextree->ctx;
	                /* This bit is tricky! For now we'll put the prob in the final link only */
                        /* pnode->logs2prob = (fsg_link_logs2prob(fsglink) >> SENSCR_SHIFT)
                           + lextree->wip + lextree->pip; */
                        pnode->logs2prob = lextree->wip + lextree->pip;
                        pnode->ci_ext = dict_first_phone(lextree->dict, dictwid);
                        pnode->ppos = 0;
                        pnode->leaf = FALSE;
                        pnode->sibling = root;  /* All root nodes linked together */
                        pnode->alloc_next = head;
                        head = pnode;
                        root = pnode;
                        ++n_lc_alloc;

                        hmm_init(lextree->ctx, &pnode->hmm, FALSE, ssid, tmatid);

                        lc_pnodelist =
                            glist_add_ptr(lc_pnodelist, (void *) pnode);
                        ssid_pnode_map[j] = pnode;
                    }
                    fsg_pnode_add_ctxt(pnode, lc);
                }
		/* Put the lc_pnodelist back into glist */
		glist->glist = lc_pnodelist;

                /* The predecessor node for the future tree is the root */
		pred = root;
            }
            else if (p != pronlen - 1) {        /* Word internal phone */
                fsg_pnode_t    *pnodeyoungest;

                ssid = dict2pid_internal(lextree->d2p, dictwid, p);
                tmatid = bin_mdef_pid2tmatid(lextree->mdef, dict_pron (lextree->dict, dictwid, p));
	        /* First check if we already have this ssid in our tree */
		pnode = pred->next.succ;
		pnodeyoungest = pnode; /* The youngest sibling */
		while (pnode && (hmm_nonmpx_ssid(&pnode->hmm) != ssid || pnode->leaf)) {
		    pnode = pnode->sibling;
		}
		if (pnode && (hmm_nonmpx_ssid(&pnode->hmm) == ssid && !pnode->leaf)) {
		    /* Found the ssid; go to next phoneme */
                    E_DEBUG(2,("Found match for %d\n", ci));
		    pred = pnode;
		    continue;
		}

		/* pnode not found, allocate it */
                pnode = (fsg_pnode_t *) ckd_calloc(1, sizeof(fsg_pnode_t));
                pnode->ctx = lextree->ctx;
                pnode->logs2prob = lextree->pip;
                pnode->ci_ext = dict_pron(lextree->dict, dictwid, p);
                pnode->ppos = p;
                pnode->leaf = FALSE;
                pnode->sibling = pnodeyoungest; /* May be NULL */
                if (p == 1) {   /* Predecessor = set of root nodes for left ctxts */
                    for (gn = lc_pnodelist; gn; gn = gnode_next(gn)) {
                        pred = (fsg_pnode_t *) gnode_ptr(gn);
                        pred->next.succ = pnode;
                    }
                }
                else {          /* Predecessor = word internal node */
                    pred->next.succ = pnode;
                }
                pnode->alloc_next = head;
                head = pnode;
                ++n_int_alloc;

                hmm_init(lextree->ctx, &pnode->hmm, FALSE, ssid, tmatid);

                pred = pnode;
            }
            else {              /* Leaf phone, handle required right contexts */
	        /* Note, leaf phones are not part of the tree */
                xwdssid_t *rssid;
                memset((void *) ssid_pnode_map, 0,
                       n_ci * sizeof(fsg_pnode_t *));
                lc = dict_pron(lextree->dict, dictwid, p-1);
                rssid = dict2pid_rssid(lextree->d2p, ci, lc);
                tmatid = bin_mdef_pid2tmatid(lextree->mdef, dict_pron (lextree->dict, dictwid, p));

                for (i = 0; rclist[i] >= 0; i++) {
                    rc = rclist[i];

                    j = rssid->cimap[rc];
                    ssid = rssid->ssid[j];
                    pnode = ssid_pnode_map[j];

                    if (!pnode) {       /* Allocate pnode for this new ssid */
                        pnode =
                            (fsg_pnode_t *) ckd_calloc(1,
                                                       sizeof
                                                       (fsg_pnode_t));
                        pnode->ctx = lextree->ctx;
			/* We are plugging the word prob here. Ugly */
                        /* pnode->logs2prob = lextree->pip; */
                        pnode->logs2prob = (fsg_link_logs2prob(fsglink) >> SENSCR_SHIFT)
                            + lextree->pip;
                        pnode->ci_ext = dict_pron(lextree->dict, dictwid, p);
                        pnode->ppos = p;
                        pnode->leaf = TRUE;
                        pnode->sibling = rc_pnodelist ?
                            (fsg_pnode_t *) gnode_ptr(rc_pnodelist) : NULL;
                        pnode->next.fsglink = fsglink;
                        pnode->alloc_next = head;
                        head = pnode;
                        ++n_rc_alloc;

                        hmm_init(lextree->ctx, &pnode->hmm, FALSE, ssid, tmatid);

                        rc_pnodelist =
                            glist_add_ptr(rc_pnodelist, (void *) pnode);
                        ssid_pnode_map[j] = pnode;
                    }
                    else {
                        assert(hmm_nonmpx_ssid(&pnode->hmm) == ssid);
                    }
                    fsg_pnode_add_ctxt(pnode, rc);
                }

                if (p == 1) {   /* Predecessor = set of root nodes for left ctxts */
                    for (gn = lc_pnodelist; gn; gn = gnode_next(gn)) {
                        pred = (fsg_pnode_t *) gnode_ptr(gn);
                        if (!pred->next.succ)
                            pred->next.succ = (fsg_pnode_t *) gnode_ptr(rc_pnodelist);
                        else {
                            /* Link to the end of the sibling chain */
                            fsg_pnode_t *succ = pred->next.succ;
                            while (succ->sibling) succ = succ->sibling;
                            succ->sibling = (fsg_pnode_t*) gnode_ptr(rc_pnodelist);
                            /* Since all entries of lc_pnodelist point
                               to the same array, sufficient to update it once */
                            break; 
                        }
                    }
                }
                else {          /* Predecessor = word internal node */
                    if (!pred->next.succ)
                        pred->next.succ = (fsg_pnode_t *) gnode_ptr(rc_pnodelist);
                    else {
                        /* Link to the end of the sibling chain */
                        fsg_pnode_t *succ = pred->next.succ;
                        while (succ->sibling) succ = succ->sibling;
                        succ->sibling = (fsg_pnode_t *) gnode_ptr(rc_pnodelist);
                    }
                }
            }
        }

        ckd_free((void *) ssid_pnode_map);
        /* glist_free(lc_pnodelist);  Nope; this gets freed outside */
        glist_free(rc_pnodelist);
    }

    E_DEBUG(2,("Allocated %d HMMs (%d lc, %d rc, %d internal)\n",
               n_lc_alloc + n_rc_alloc + n_int_alloc,
               n_lc_alloc, n_rc_alloc, n_int_alloc));
    *alloc_head = head;

    return root;
}


static fsg_pnode_t *
fsg_psubtree_init(fsg_lextree_t *lextree,
                  fsg_model_t * fsg, int32 from_state,
                  fsg_pnode_t ** alloc_head)
{
    int32 dst;
    gnode_t *gn;
    fsg_link_t *fsglink;
    fsg_pnode_t *root;
    int32 n_ci, n_arc;
    fsg_glist_linklist_t *glist = NULL;

    root = NULL;
    assert(*alloc_head == NULL);

    n_ci = bin_mdef_n_ciphone(lextree->mdef);
    if (n_ci > (FSG_PNODE_CTXT_BVSZ * 32)) {
        E_FATAL
            ("#phones > %d; increase FSG_PNODE_CTXT_BVSZ and recompile\n",
             FSG_PNODE_CTXT_BVSZ * 32);
    }
    n_arc = 0;
    for (dst = 0; dst < fsg_model_n_state(fsg); dst++) {
        /* Add all links from from_state to dst */
        for (gn = fsg_model_trans(fsg, from_state, dst); gn;
             gn = gnode_next(gn)) {
            /* Add word emitted by this transition (fsglink) to lextree */
            fsglink = (fsg_link_t *) gnode_ptr(gn);

            assert(fsg_link_wid(fsglink) >= 0);     /* Cannot be a null trans */

            E_DEBUG(2,("Building lextree for arc from %d to %d: %s\n",
                       from_state, dst, fsg_model_word_str(fsg, fsg_link_wid(fsglink))));
            root = psubtree_add_trans(lextree, root, &glist, fsglink,
                                      lextree->lc[from_state],
                                      lextree->rc[dst],
                                      alloc_head);
            ++n_arc;
        }
    }
    E_DEBUG(2,("State %d has %d outgoing arcs\n", from_state, n_arc));

    fsg_glist_linklist_free(glist);

    return root;
}


static void
fsg_psubtree_free(fsg_pnode_t * head)
{
    fsg_pnode_t *next;

    while (head) {
        next = head->alloc_next;
        hmm_deinit(&head->hmm);
        ckd_free(head);
        head = next;
    }
}

void fsg_psubtree_dump_node(fsg_lextree_t *tree, fsg_pnode_t *node, FILE *fp)
{    
    int32 i;
    fsg_link_t *tl;

    /* Indentation */
    for (i = 0; i <= node->ppos; i++)
        fprintf(fp, "  ");

    fprintf(fp, "%p.@", node);    /* Pointer used as node
                    		   * ID */
    fprintf(fp, " %5d.SS", hmm_nonmpx_ssid(&node->hmm));
    fprintf(fp, " %10d.LP", node->logs2prob);
    fprintf(fp, " %p.SIB", node->sibling);
    fprintf(fp, " %s.%d", bin_mdef_ciphone_str(tree->mdef, node->ci_ext), node->ppos);
    if ((node->ppos == 0) || node->leaf) {
        fprintf(fp, " [");
        for (i = 0; i < FSG_PNODE_CTXT_BVSZ; i++)
            fprintf(fp, "%08x", node->ctxt.bv[i]);
        fprintf(fp, "]");
    }
    if (node->leaf) {
        tl = node->next.fsglink;
        fprintf(fp, " {%s[%d->%d](%d)}",
                fsg_model_word_str(tree->fsg, tl->wid),
                tl->from_state, tl->to_state, tl->logs2prob);
    } else {
        fprintf(fp, " %p.NXT", node->next.succ);
    }
    fprintf(fp, "\n");

    return;
}

void 
fsg_psubtree_dump(fsg_lextree_t *tree, fsg_pnode_t *root, FILE * fp)
{
    fsg_pnode_t *succ;

    if (root == NULL) return;
    if (root->ppos == 0) {
        while(root->sibling && root->sibling->next.succ == root->next.succ) {
            fsg_psubtree_dump_node(tree, root, fp);
            root = root->sibling;
        }
        fflush(fp);
    }
    
    fsg_psubtree_dump_node(tree, root, fp);

    if (root->leaf) {
        if (root->ppos == 0 && root->sibling) { // For single-phone words
            fsg_psubtree_dump(tree, root->sibling,fp);
        }
        return;
    }

    succ = root->next.succ;
    while(succ) {
        fsg_psubtree_dump(tree, succ,fp);
        succ = succ->sibling;
    }

    if (root->ppos == 0) {
        fsg_psubtree_dump(tree, root->sibling,fp);
        fflush(fp);
    }

    return;
}

void
fsg_psubtree_pnode_deactivate(fsg_pnode_t * pnode)
{
    hmm_clear(&pnode->hmm);
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * fsg_lextree.h -- The collection of all the lextrees for the entire FSM.
 * 
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 2004 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * 
 * $Log: fsg_lextree.h,v $
 * Revision 1.1.1.1  2006/05/23 18:45:02  dhuggins
 * re-importation
 *
 * Revision 1.1  2004/07/16 00:57:12  egouvea
 * Added Ravi's implementation of FSG support.
 *
 * Revision 1.3  2004/06/23 20:32:16  rkm
 * *** empty log message ***
 *
 * Revision 1.2  2004/05/27 14:22:57  rkm
 * FSG cross-word triphones completed (but for single-phone words)
 *
 * Revision 1.1.1.1  2004/03/01 14:30:31  rkm
 *
 *
 * Revision 1.1  2004/02/23 15:53:45  rkm
 * Renamed from fst to fsg
 *
 * Revision 1.2  2004/02/19 21:16:54  rkm
 * Added fsg_search.{c,h}
 *
 * Revision 1.1  2004/02/18 15:02:34  rkm
 * Added fsg_lextree.{c,h}
 *
 * 
 * 18-Feb-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon
 * 		Started.
 */


#ifndef __S2_FSG_LEXTREE_H__
#define __S2_FSG_LEXTREE_H__

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/fsg_model.h>

/* Local headers. */
#include "hmm.h"
#include "dict.h"
#include "dict2pid.h"

/*
 * **HACK-ALERT**!!  Compile-time constant determining the size of the
 * bitvector fsg_pnode_t.fsg_pnode_ctxt_t.bv.  (See below.)
 * But it makes memory allocation simpler and more efficient.
 */
#define FSG_PNODE_CTXT_BVSZ	2

typedef struct {
    uint32 bv[FSG_PNODE_CTXT_BVSZ];
} fsg_pnode_ctxt_t;


/*
 * All transitions (words) out of any given FSG state represented are by a
 * phonetic prefix lextree (except for epsilon or null transitions; they
 * are not part of the lextree).  Lextree leaf nodes represent individual
 * FSG transitions, so no sharing is allowed at the leaf nodes.  The FSG
 * transition probs are distributed along the lextree: the prob at a node
 * is the max of the probs of all leaf nodes (and, hence, FSG transitions)
 * reachable from that node.
 * 
 * To conserve memory, the underlying HMMs with state-level information are
 * allocated only as needed.  Root and leaf nodes must also account for all
 * the possible phonetic contexts, with an independent HMM for each distinct
 * context.
 */
typedef struct fsg_pnode_s {
    /*
     * If this is not a leaf node, the first successor (child) node.  Otherwise
     * the parent FSG transition for which this is the leaf node (for figuring
     * the FSG destination state, and word emitted by the transition).  A node
     * may have several children.  The succ ptr gives just the first; the rest
     * are linked via the sibling ptr below.
     */
    union {
        struct fsg_pnode_s *succ;
        fsg_link_t *fsglink;
    } next;
  
    /*
     * For simplicity of memory management (i.e., freeing the pnodes), all
     * pnodes allocated for all transitions out of a state are maintained in a
     * linear linked list through the alloc_next pointer.
     */
    struct fsg_pnode_s *alloc_next;
  
    /*
     * The next node that is also a child of the parent of this node; NULL if
     * none.
     */
    struct fsg_pnode_s *sibling;

    /*
     * The transition (log) probability to be incurred upon transitioning to
     * this node.  (Transition probabilities are really associated with the
     * transitions.  But a lextree node has exactly one incoming transition.
     * Hence, the prob can be associated with the node.)
     * This is a logs2(prob) value, and includes the language weight.
     */
    int32 logs2prob;
  
    /*
     * The root and leaf positions associated with any transition have to deal
     * with multiple phonetic contexts.  However, different contexts may result
     * in the same SSID (senone-seq ID), and can share a single pnode with that
     * SSID.  But the pnode should track the set of context CI phones that share
     * it.  Hence the fsg_pnode_ctxt_t bit-vector set-representation.  (For
     * simplicity of implementation, its size is a compile-time constant for
     * now.)  Single phone words would need a 2-D array of context, but that's
     * too expensive.  For now, they simply use SIL as right context, so only
     * the left context is properly modelled.
     * (For word-internal phones, this field is unused, of course.)
     */
    fsg_pnode_ctxt_t ctxt;
  
    uint16 ci_ext;	/* This node's CIphone as viewed externally (context) */
    uint8 ppos;	/* Phoneme position in pronunciation */
    uint8 leaf;	/* Whether this is a leaf node */
  
    /* HMM-state-level stuff here */
    hmm_context_t *ctx;
    hmm_t hmm;
} fsg_pnode_t;

/* Access macros */
#define fsg_pnode_leaf(p)	((p)->leaf)
#define fsg_pnode_logs2prob(p)	((p)->logs2prob)
#define fsg_pnode_succ(p)	((p)->next.succ)
#define fsg_pnode_fsglink(p)	((p)->next.fsglink)
#define fsg_pnode_sibling(p)	((p)->sibling)
#define fsg_pnode_hmmptr(p)	(&((p)->hmm))
#define fsg_pnode_ci_ext(p)	((p)->ci_ext)
#define fsg_pnode_ppos(p)	((p)->ppos)
#define fsg_pnode_leaf(p)	((p)->leaf)
#define fsg_pnode_ctxt(p)	((p)->ctxt)

#define fsg_pnode_add_ctxt(p,c)	((p)->ctxt.bv[(c)>>5] |= (1 << ((c)&0x001f)))

/*
 * The following is macroized because its called very frequently
 * ::: uint32 fsg_pnode_ctxt_sub (fsg_pnode_ctxt_t *src, fsg_pnode_ctxt_t *sub);
 */
/*
 * Subtract bitvector sub from bitvector src (src updated with the result).
 * Return 0 if result is all 0, non-zero otherwise.
 */

#if (FSG_PNODE_CTXT_BVSZ == 1)
    #define FSG_PNODE_CTXT_SUB(src,sub) \
    ((src)->bv[0] = (~((sub)->bv[0]) & (src)->bv[0]))
#elif (FSG_PNODE_CTXT_BVSZ == 2)
    #define FSG_PNODE_CTXT_SUB(src,sub) \
    (((src)->bv[0] = (~((sub)->bv[0]) & (src)->bv[0])) | \
     ((src)->bv[1] = (~((sub)->bv[1]) & (src)->bv[1])))
#elif (FSG_PNODE_CTXT_BVSZ == 4)
    #define FSG_PNODE_CTXT_SUB(src,sub) \
    (((src)->bv[0] = (~((sub)->bv[0]) & (src)->bv[0]))  | \
     ((src)->bv[1] = (~((sub)->bv[1]) & (src)->bv[1]))  | \
     ((src)->bv[2] = (~((sub)->bv[2]) & (src)->bv[2]))  | \
     ((src)->bv[3] = (~((sub)->bv[3]) & (src)->bv[3])))
#else
    #define FSG_PNODE_CTXT_SUB(src,sub) fsg_pnode_ctxt_sub_generic((src),(sub))
#endif

/**
 * Collection of lextrees for an FSG.
 */
typedef struct fsg_lextree_s {
    fsg_model_t *fsg;	/**< The fsg for which this lextree is built. */
    hmm_context_t *ctx; /**< HMM context structure. */
    dict_t *dict;     /**< Pronunciation dictionary for this FSG. */
    dict2pid_t *d2p;    /**< Context-dependent phone mappings for this FSG. */
    bin_mdef_t *mdef;   /**< Model definition (triphone mappings). */

    /*
     * Left and right CIphone sets for each state.
     * Left context CIphones for a state S: If word W transitions into S, W's
     * final CIphone is in S's {lc}.  Words transitioning out of S must consider
     * these left context CIphones.
     * Similarly, right contexts for state S: If word W transitions out of S,
     * W's first CIphone is in S's {rc}.  Words transitioning into S must consider
     * these right contexts.
     * 
     * NOTE: Words may transition into and out of S INDIRECTLY, with intermediate
     *   null transitions.
     * NOTE: Single-phone words are difficult; only SILENCE right context is
     *   modelled for them.
     * NOTE: Non-silence filler phones aren't included in these sets.  Filler
     *   words don't use context, and present the SILENCE phone as context to
     *   adjacent words.
     */
    int16 **lc;         /**< Left context triphone mappings for FSG. */
    int16 **rc;         /**< Right context triphone mappings for FSG. */

    fsg_pnode_t **root;	/* root[s] = lextree representing all transitions
			   out of state s.  Note that the "tree" for each
			   state is actually a collection of trees, linked
			   via fsg_pnode_t.sibling (root[s]->sibling) */
    fsg_pnode_t **alloc_head;	/* alloc_head[s] = head of linear list of all
				   pnodes allocated for state s */
    int32 n_pnode;	/* #HMM nodes in search structure */
    int32 wip;
    int32 pip;
} fsg_lextree_t;

/* Access macros */
#define fsg_lextree_root(lt,s)	((lt)->root[s])
#define fsg_lextree_n_pnode(lt)	((lt)->n_pnode)

/**
 * Create, initialize, and return a new phonetic lextree for the given FSG.
 */
fsg_lextree_t *fsg_lextree_init(fsg_model_t *fsg, dict_t *dict,
                                dict2pid_t *d2p,
				bin_mdef_t *mdef, hmm_context_t *ctx,
				int32 wip, int32 pip);

/**
 * Free lextrees for an FSG.
 */
void fsg_lextree_free(fsg_lextree_t *fsg);

/**
 * Print an FSG lextree to a file for debugging.
 */
void fsg_lextree_dump(fsg_lextree_t *fsg, FILE *fh);

/**
 * Mark the given pnode as inactive (for search).
 */
void fsg_psubtree_pnode_deactivate(fsg_pnode_t *pnode);

/**
 *  Set all flags on in the given context bitvector.
 */
void fsg_pnode_add_all_ctxt(fsg_pnode_ctxt_t *ctxt);

/**
 *  Generic variant for arbitrary size
 */
uint32 fsg_pnode_ctxt_sub_generic(fsg_pnode_ctxt_t *src, fsg_pnode_ctxt_t *sub);

#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * fsg_search.c -- Search structures for FSM decoding.
 * 
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 2004 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 *
 * 18-Feb-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon
 * 		Started.
 */

/* System headers. */
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/cmd_ln.h>

/* Local headers. */
#include "pocketsphinx_internal.h"
#include "ps_lattice_internal.h"
#include "fsg_search_internal.h"
#include "fsg_history.h"
#include "fsg_lextree.h"

/* Turn this on for detailed debugging dump */
#define __FSG_DBG__		0
#define __FSG_DBG_CHAN__	0

static ps_seg_t *fsg_search_seg_iter(ps_search_t *search, int32 *out_score);
static ps_lattice_t *fsg_search_lattice(ps_search_t *search);
static int fsg_search_prob(ps_search_t *search);

static ps_searchfuncs_t fsg_funcs = {
    /* name: */   "fsg",
    /* start: */  fsg_search_start,
    /* step: */   fsg_search_step,
    /* finish: */ fsg_search_finish,
    /* reinit: */ fsg_search_reinit,
    /* free: */   fsg_search_free,
    /* lattice: */  fsg_search_lattice,
    /* hyp: */      fsg_search_hyp,
    /* prob: */     fsg_search_prob,
    /* seg_iter: */ fsg_search_seg_iter,
};

ps_search_t *
fsg_search_init(cmd_ln_t *config,
                acmod_t *acmod,
                dict_t *dict,
                dict2pid_t *d2p)
{
    fsg_search_t *fsgs;
    char const *path;

    fsgs = ckd_calloc(1, sizeof(*fsgs));
    ps_search_init(ps_search_base(fsgs), &fsg_funcs, config, acmod, dict, d2p);

    /* Initialize HMM context. */
    fsgs->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
                                    acmod->tmat->tp, NULL, acmod->mdef->sseq);
    if (fsgs->hmmctx == NULL) {
        ps_search_free(ps_search_base(fsgs));
        return NULL;
    }

    /* Intialize the search history object */
    fsgs->history = fsg_history_init(NULL, dict);
    fsgs->frame = -1;

    /* Initialize FSG table. */
    fsgs->fsgs = hash_table_new(5, HASH_CASE_YES);

    /* Get search pruning parameters */
    fsgs->beam_factor = 1.0f;
    fsgs->beam = fsgs->beam_orig
        = (int32) logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-beam"))
        >> SENSCR_SHIFT;
    fsgs->pbeam = fsgs->pbeam_orig
        = (int32) logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-pbeam"))
        >> SENSCR_SHIFT;
    fsgs->wbeam = fsgs->wbeam_orig
        = (int32) logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-wbeam"))
        >> SENSCR_SHIFT;

    /* LM related weights/penalties */
    fsgs->lw = cmd_ln_float32_r(config, "-lw");
    fsgs->pip = (int32) (logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-pip"))
                           * fsgs->lw)
        >> SENSCR_SHIFT;
    fsgs->wip = (int32) (logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-wip"))
                           * fsgs->lw)
        >> SENSCR_SHIFT;

    /* Best path search (and confidence annotation)? */
    if (cmd_ln_boolean_r(config, "-bestpath"))
        fsgs->bestpath = TRUE;

    /* Acoustic score scale for posterior probabilities. */
    fsgs->ascale = 1.0 / cmd_ln_float32_r(config, "-ascale");

    E_INFO("FSG(beam: %d, pbeam: %d, wbeam: %d; wip: %d, pip: %d)\n",
           fsgs->beam_orig, fsgs->pbeam_orig, fsgs->wbeam_orig,
           fsgs->wip, fsgs->pip);

    /* Load an FSG if one was specified in config */
    if ((path = cmd_ln_str_r(config, "-fsg"))) {
        fsg_model_t *fsg;

        if ((fsg = fsg_model_readfile(path, acmod->lmath, fsgs->lw)) == NULL)
            goto error_out;
        if (fsg_set_add(fsgs, fsg_model_name(fsg), fsg) != fsg) {
            fsg_model_free(fsg);
            goto error_out;
        }
        if (fsg_set_select(fsgs, fsg_model_name(fsg)) == NULL)
            goto error_out;
        if (fsg_search_reinit(ps_search_base(fsgs),
                              ps_search_dict(fsgs),
                              ps_search_dict2pid(fsgs)) < 0)
            goto error_out;
    }
    /* Or load a JSGF grammar */
    else if ((path = cmd_ln_str_r(config, "-jsgf"))) {
        fsg_model_t *fsg;
        jsgf_rule_t *rule;
        char const *toprule;

        if ((fsgs->jsgf = jsgf_parse_file(path, NULL)) == NULL)
            goto error_out;

        rule = NULL;
        /* Take the -toprule if specified. */
        if ((toprule = cmd_ln_str_r(config, "-toprule"))) {
            char *anglerule;
            anglerule = string_join("<", toprule, ">", NULL);
            rule = jsgf_get_rule(fsgs->jsgf, anglerule);
            ckd_free(anglerule);
            if (rule == NULL) {
                E_ERROR("Start rule %s not found\n", toprule);
                goto error_out;
            }
        }
        /* Otherwise, take the first public rule. */
        else {
            jsgf_rule_iter_t *itor;

            for (itor = jsgf_rule_iter(fsgs->jsgf); itor;
                 itor = jsgf_rule_iter_next(itor)) {
                rule = jsgf_rule_iter_rule(itor);
                if (jsgf_rule_public(rule)) {
            	    jsgf_rule_iter_free(itor);
                    break;
                }
            }
            if (rule == NULL) {
                E_ERROR("No public rules found in %s\n", path);
                goto error_out;
            }
        }
        fsg = jsgf_build_fsg(fsgs->jsgf, rule, acmod->lmath, fsgs->lw);
        if (fsg_set_add(fsgs, fsg_model_name(fsg), fsg) != fsg) {
            fsg_model_free(fsg);
            goto error_out;
        }
        if (fsg_set_select(fsgs, fsg_model_name(fsg)) == NULL)
            goto error_out;
        if (fsg_search_reinit(ps_search_base(fsgs),
                              ps_search_dict(fsgs),
                              ps_search_dict2pid(fsgs)) < 0)
            goto error_out;
    }

    return ps_search_base(fsgs);

error_out:
    fsg_search_free(ps_search_base(fsgs));
    return NULL;
}

void
fsg_search_free(ps_search_t *search)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;
    hash_iter_t *itor;

    ps_search_deinit(search);
    if (fsgs->jsgf)
        jsgf_grammar_free(fsgs->jsgf);
    fsg_lextree_free(fsgs->lextree);
    if (fsgs->history) {
        fsg_history_reset(fsgs->history);
        fsg_history_set_fsg(fsgs->history, NULL, NULL);
        fsg_history_free(fsgs->history);
    }
    if (fsgs->fsgs) {
        for (itor = hash_table_iter(fsgs->fsgs);
             itor; itor = hash_table_iter_next(itor)) {
            fsg_model_t *fsg = (fsg_model_t *) hash_entry_val(itor->ent);
            fsg_model_free(fsg);
        }
        hash_table_free(fsgs->fsgs);
    }
    hmm_context_free(fsgs->hmmctx);
    ckd_free(fsgs);
}

int
fsg_search_reinit(ps_search_t *search, dict_t *dict, dict2pid_t *d2p)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;

    /* Free the old lextree */
    if (fsgs->lextree)
        fsg_lextree_free(fsgs->lextree);

    /* Free old dict2pid, dict */
    ps_search_base_reinit(search, dict, d2p);
    
    /* Nothing to update */
    if (fsgs->fsg == NULL)
	return 0;

    /* Update the number of words (not used by this module though). */
    search->n_words = dict_size(dict);

    /* Allocate new lextree for the given FSG */
    fsgs->lextree = fsg_lextree_init(fsgs->fsg, dict, d2p,
                                     ps_search_acmod(fsgs)->mdef,
                                     fsgs->hmmctx, fsgs->wip, fsgs->pip);

    /* Inform the history module of the new fsg */
    fsg_history_set_fsg(fsgs->history, fsgs->fsg, dict);

    return 0;
}


static int
fsg_search_add_silences(fsg_search_t *fsgs, fsg_model_t *fsg)
{
    dict_t *dict;
    int32 wid;
    int n_sil;

    dict = ps_search_dict(fsgs);
    /*
     * NOTE: Unlike N-Gram search, we do not use explicit start and
     * end symbols.  This is because the start and end nodes are
     * defined in the grammar.  We do add silence/filler self-loops to
     * all states in order to allow for silence between words and at
     * the beginning and end of utterances.
     *
     * This has some implications for word graph generation, namely,
     * that there can (and usually will) be multiple start and end
     * states in the word graph.  We therefore do add explicit start
     * and end nodes to the graph.
     */
    /* Add silence self-loops to all states. */
    fsg_model_add_silence(fsg, "<sil>", -1,
                          cmd_ln_float32_r(ps_search_config(fsgs), "-silprob"));
    n_sil = 0;
    /* Add self-loops for all other fillers. */
    for (wid = dict_filler_start(dict); wid < dict_filler_end(dict); ++wid) {
        char const *word = dict_wordstr(dict, wid);
        if (wid == dict_startwid(dict) || wid == dict_finishwid(dict))
            continue;
        fsg_model_add_silence(fsg, word, -1,
                              cmd_ln_float32_r(ps_search_config(fsgs), "-fillprob"));
        ++n_sil;
    }

    return n_sil;
}

/* Scans the dictionary and check if all words are present. */
static int
fsg_search_check_dict(fsg_search_t *fsgs, fsg_model_t *fsg)
{
    dict_t *dict;
    int i;

    dict = ps_search_dict(fsgs);
    for (i = 0; i < fsg_model_n_word(fsg); ++i) {
        char const *word;
        int32 wid;

        word = fsg_model_word_str(fsg, i);
        wid = dict_wordid(dict, word);
        if (wid == BAD_S3WID) {
    	    E_ERROR("The word '%s' is missing in the dictionary\n", word);
    	    return FALSE;
    	}
    }

    return TRUE;
}

static int
fsg_search_add_altpron(fsg_search_t *fsgs, fsg_model_t *fsg)
{
    dict_t *dict;
    int n_alt, n_word;
    int i;

    dict = ps_search_dict(fsgs);
    /* Scan FSG's vocabulary for words that have alternate pronunciations. */
    n_alt = 0;
    n_word = fsg_model_n_word(fsg);
    for (i = 0; i < n_word; ++i) {
        char const *word;
        int32 wid;

        word = fsg_model_word_str(fsg, i);
        wid = dict_wordid(dict, word);
        if (wid != BAD_S3WID) {
            while ((wid = dict_nextalt(dict, wid)) != BAD_S3WID) {
	        n_alt += fsg_model_add_alt(fsg, word, dict_wordstr(dict, wid));
    	    }
    	}
    }

    E_INFO("Added %d alternate word transitions\n", n_alt);
    return n_alt;
}

fsg_model_t *
fsg_set_get_fsg(fsg_search_t *fsgs, const char *name)
{
    void *val;

    if (hash_table_lookup(fsgs->fsgs, name, &val) < 0)
        return NULL;
    return (fsg_model_t *)val;
}

fsg_model_t *
fsg_set_add(fsg_search_t *fsgs, char const *name, fsg_model_t *fsg)
{
    if (name == NULL)
        name = fsg_model_name(fsg);

    if (!fsg_search_check_dict(fsgs, fsg))
	return NULL;

    /* Add silence transitions and alternate words. */
    if (cmd_ln_boolean_r(ps_search_config(fsgs), "-fsgusefiller")
        && !fsg_model_has_sil(fsg))
        fsg_search_add_silences(fsgs, fsg);
    if (cmd_ln_boolean_r(ps_search_config(fsgs), "-fsgusealtpron")
        && !fsg_model_has_alt(fsg))
        fsg_search_add_altpron(fsgs, fsg);

    return (fsg_model_t *)hash_table_enter(fsgs->fsgs, name, fsg);
}


fsg_model_t *
fsg_set_remove_byname(fsg_search_t *fsgs, char const *key)
{
    fsg_model_t *oldfsg;
    void *val;

    /* Look for the matching FSG. */
    if (hash_table_lookup(fsgs->fsgs, key, &val) < 0) {
        E_ERROR("FSG `%s' to be deleted not found\n", key);
        return NULL;
    }
    oldfsg = val;

    /* Remove it from the FSG table. */
    hash_table_delete(fsgs->fsgs, key);
    /* If this was the currently active FSG, also delete other stuff */
    if (fsgs->fsg == oldfsg) {
        fsg_lextree_free(fsgs->lextree);
        fsgs->lextree = NULL;
        fsg_history_set_fsg(fsgs->history, NULL, NULL);
        fsgs->fsg = NULL;
    }
    return oldfsg;
}


fsg_model_t *
fsg_set_remove(fsg_search_t *fsgs, fsg_model_t *fsg)
{
    char const *key;
    hash_iter_t *itor;

    key = NULL;
    for (itor = hash_table_iter(fsgs->fsgs);
         itor; itor = hash_table_iter_next(itor)) {
        fsg_model_t *oldfsg;

        oldfsg = (fsg_model_t *) hash_entry_val(itor->ent);
        if (oldfsg == fsg) {
            key = hash_entry_key(itor->ent);
            hash_table_iter_free(itor);
            break;
        }
    }
    if (key == NULL) {
        E_WARN("FSG '%s' to be deleted not found\n", fsg_model_name(fsg));
        return NULL;
    }
    else
        return fsg_set_remove_byname(fsgs, key);
}


fsg_model_t *
fsg_set_select(fsg_search_t *fsgs, const char *name)
{
    fsg_model_t *fsg;

    fsg = fsg_set_get_fsg(fsgs, name);
    if (fsg == NULL) {
        E_ERROR("FSG '%s' not known; cannot make it current\n", name);
        return NULL;
    }
    fsgs->fsg = fsg;
    return fsg;
}

fsg_set_iter_t *
fsg_set_iter(fsg_set_t *fsgs)
{
    return hash_table_iter(fsgs->fsgs);
}

fsg_set_iter_t *
fsg_set_iter_next(fsg_set_iter_t *itor)
{
    return hash_table_iter_next(itor);
}

fsg_model_t *
fsg_set_iter_fsg(fsg_set_iter_t *itor)
{
    return ((fsg_model_t *)itor->ent->val);
}

void
fsg_set_iter_free(fsg_set_iter_t *itor)
{
    hash_table_iter_free(itor);
}

static void
fsg_search_sen_active(fsg_search_t *fsgs)
{
    gnode_t *gn;
    fsg_pnode_t *pnode;
    hmm_t *hmm;

    acmod_clear_active(ps_search_acmod(fsgs));

    for (gn = fsgs->pnode_active; gn; gn = gnode_next(gn)) {
        pnode = (fsg_pnode_t *) gnode_ptr(gn);
        hmm = fsg_pnode_hmmptr(pnode);
        assert(hmm_frame(hmm) == fsgs->frame);
        acmod_activate_hmm(ps_search_acmod(fsgs), hmm);
    }
}


/*
 * Evaluate all the active HMMs.
 * (Executed once per frame.)
 */
static void
fsg_search_hmm_eval(fsg_search_t *fsgs)
{
    gnode_t *gn;
    fsg_pnode_t *pnode;
    hmm_t *hmm;
    int32 bestscore;
    int32 n, maxhmmpf;

    bestscore = WORST_SCORE;

    if (!fsgs->pnode_active) {
        E_ERROR("Frame %d: No active HMM!!\n", fsgs->frame);
        return;
    }

    for (n = 0, gn = fsgs->pnode_active; gn; gn = gnode_next(gn), n++) {
        int32 score;

        pnode = (fsg_pnode_t *) gnode_ptr(gn);
        hmm = fsg_pnode_hmmptr(pnode);
        assert(hmm_frame(hmm) == fsgs->frame);

#if __FSG_DBG__
        E_INFO("pnode(%08x) active @frm %5d\n", (int32) pnode,
               fsgs->frame);
        hmm_dump(hmm, stdout);
#endif
        score = hmm_vit_eval(hmm);
#if __FSG_DBG_CHAN__
        E_INFO("pnode(%08x) after eval @frm %5d\n",
               (int32) pnode, fsgs->frame);
        hmm_dump(hmm, stdout);
#endif

        if (score BETTER_THAN bestscore)
            bestscore = score;
    }

#if __FSG_DBG__
    E_INFO("[%5d] %6d HMM; bestscr: %11d\n", fsgs->frame, n, bestscore);
#endif
    fsgs->n_hmm_eval += n;

    /* Adjust beams if #active HMMs larger than absolute threshold */
    maxhmmpf = cmd_ln_int32_r(ps_search_config(fsgs), "-maxhmmpf");
    if (maxhmmpf != -1 && n > maxhmmpf) {
        /*
         * Too many HMMs active; reduce the beam factor applied to the default
         * beams, but not if the factor is already at a floor (0.1).
         */
        if (fsgs->beam_factor > 0.1) {        /* Hack!!  Hardwired constant 0.1 */
            fsgs->beam_factor *= 0.9f;        /* Hack!!  Hardwired constant 0.9 */
            fsgs->beam =
                (int32) (fsgs->beam_orig * fsgs->beam_factor);
            fsgs->pbeam =
                (int32) (fsgs->pbeam_orig * fsgs->beam_factor);
            fsgs->wbeam =
                (int32) (fsgs->wbeam_orig * fsgs->beam_factor);
        }
    }
    else {
        fsgs->beam_factor = 1.0f;
        fsgs->beam = fsgs->beam_orig;
        fsgs->pbeam = fsgs->pbeam_orig;
        fsgs->wbeam = fsgs->wbeam_orig;
    }

    if (n > fsg_lextree_n_pnode(fsgs->lextree))
        E_FATAL("PANIC! Frame %d: #HMM evaluated(%d) > #PNodes(%d)\n",
                fsgs->frame, n, fsg_lextree_n_pnode(fsgs->lextree));

    fsgs->bestscore = bestscore;
}


static void
fsg_search_pnode_trans(fsg_search_t *fsgs, fsg_pnode_t * pnode)
{
    fsg_pnode_t *child;
    hmm_t *hmm;
    int32 newscore, thresh, nf;

    assert(pnode);
    assert(!fsg_pnode_leaf(pnode));

    nf = fsgs->frame + 1;
    thresh = fsgs->bestscore + fsgs->beam;

    hmm = fsg_pnode_hmmptr(pnode);

    for (child = fsg_pnode_succ(pnode);
         child; child = fsg_pnode_sibling(child)) {
        newscore = hmm_out_score(hmm) + child->logs2prob;

        if ((newscore BETTER_THAN thresh)
            && (newscore BETTER_THAN hmm_in_score(&child->hmm))) {
            /* Incoming score > pruning threshold and > target's existing score */
            if (hmm_frame(&child->hmm) < nf) {
                /* Child node not yet activated; do so */
                fsgs->pnode_active_next =
                    glist_add_ptr(fsgs->pnode_active_next,
                                  (void *) child);
            }

            hmm_enter(&child->hmm, newscore, hmm_out_history(hmm), nf);
        }
    }
}


static void
fsg_search_pnode_exit(fsg_search_t *fsgs, fsg_pnode_t * pnode)
{
    hmm_t *hmm;
    fsg_link_t *fl;
    int32 wid;
    fsg_pnode_ctxt_t ctxt;

    assert(pnode);
    assert(fsg_pnode_leaf(pnode));

    hmm = fsg_pnode_hmmptr(pnode);
    fl = fsg_pnode_fsglink(pnode);
    assert(fl);

    wid = fsg_link_wid(fl);
    assert(wid >= 0);

#if __FSG_DBG__
    E_INFO("[%5d] Exit(%08x) %10d(score) %5d(pred)\n",
           fsgs->frame, (int32) pnode,
           hmm_out_score(hmm), hmm_out_history(hmm));
#endif

    /*
     * Check if this is filler or single phone word; these do not model right
     * context (i.e., the exit score applies to all right contexts).
     */
    if (fsg_model_is_filler(fsgs->fsg, wid)
        /* FIXME: This might be slow due to repeated calls to dict_to_id(). */
        || (dict_is_single_phone(ps_search_dict(fsgs),
                                   dict_wordid(ps_search_dict(fsgs),
                                               fsg_model_word_str(fsgs->fsg, wid))))) {
        /* Create a dummy context structure that applies to all right contexts */
        fsg_pnode_add_all_ctxt(&ctxt);

        /* Create history table entry for this word exit */
        fsg_history_entry_add(fsgs->history,
                              fl,
                              fsgs->frame,
                              hmm_out_score(hmm),
                              hmm_out_history(hmm),
                              pnode->ci_ext, ctxt);

    }
    else {
        /* Create history table entry for this word exit */
        fsg_history_entry_add(fsgs->history,
                              fl,
                              fsgs->frame,
                              hmm_out_score(hmm),
                              hmm_out_history(hmm),
                              pnode->ci_ext, pnode->ctxt);
    }
}


/*
 * (Beam) prune the just evaluated HMMs, determine which ones remain
 * active, which ones transition to successors, which ones exit and
 * terminate in their respective destination FSM states.
 * (Executed once per frame.)
 */
static void
fsg_search_hmm_prune_prop(fsg_search_t *fsgs)
{
    gnode_t *gn;
    fsg_pnode_t *pnode;
    hmm_t *hmm;
    int32 thresh, word_thresh, phone_thresh;

    assert(fsgs->pnode_active_next == NULL);

    thresh = fsgs->bestscore + fsgs->beam;
    phone_thresh = fsgs->bestscore + fsgs->pbeam;
    word_thresh = fsgs->bestscore + fsgs->wbeam;

    for (gn = fsgs->pnode_active; gn; gn = gnode_next(gn)) {
        pnode = (fsg_pnode_t *) gnode_ptr(gn);
        hmm = fsg_pnode_hmmptr(pnode);

        if (hmm_bestscore(hmm) >= thresh) {
            /* Keep this HMM active in the next frame */
            if (hmm_frame(hmm) == fsgs->frame) {
                hmm_frame(hmm) = fsgs->frame + 1;
                fsgs->pnode_active_next =
                    glist_add_ptr(fsgs->pnode_active_next,
                                  (void *) pnode);
            }
            else {
                assert(hmm_frame(hmm) == fsgs->frame + 1);
            }

            if (!fsg_pnode_leaf(pnode)) {
                if (hmm_out_score(hmm) >= phone_thresh) {
                    /* Transition out of this phone into its children */
                    fsg_search_pnode_trans(fsgs, pnode);
                }
            }
            else {
                if (hmm_out_score(hmm) >= word_thresh) {
                    /* Transition out of leaf node into destination FSG state */
                    fsg_search_pnode_exit(fsgs, pnode);
                }
            }
        }
    }
}


/*
 * Propagate newly created history entries through null transitions.
 */
static void
fsg_search_null_prop(fsg_search_t *fsgs)
{
    int32 bpidx, n_entries, thresh, newscore;
    fsg_hist_entry_t *hist_entry;
    fsg_link_t *l;
    int32 s;
    fsg_model_t *fsg;

    fsg = fsgs->fsg;
    thresh = fsgs->bestscore + fsgs->wbeam; /* Which beam really?? */

    n_entries = fsg_history_n_entries(fsgs->history);

    for (bpidx = fsgs->bpidx_start; bpidx < n_entries; bpidx++) {
        fsg_arciter_t *itor;
        hist_entry = fsg_history_entry_get(fsgs->history, bpidx);

        l = fsg_hist_entry_fsglink(hist_entry);

        /* Destination FSG state for history entry */
        s = l ? fsg_link_to_state(l) : fsg_model_start_state(fsg);

        /*
         * Check null transitions from d to all other states.  (Only need to
         * propagate one step, since FSG contains transitive closure of null
         * transitions.)
         */
        /* Add all links from from_state to dst */
        for (itor = fsg_model_arcs(fsg, s); itor;
             itor = fsg_arciter_next(itor)) {
            fsg_link_t *l = fsg_arciter_get(itor);

            /* FIXME: Need to deal with tag transitions somehow. */
            if (fsg_link_wid(l) != -1)
                continue;
            newscore =
                fsg_hist_entry_score(hist_entry) +
                (fsg_link_logs2prob(l) >> SENSCR_SHIFT);

            if (newscore >= thresh) {
                fsg_history_entry_add(fsgs->history, l,
                                      fsg_hist_entry_frame(hist_entry),
                                      newscore,
                                      bpidx,
                                      fsg_hist_entry_lc(hist_entry),
                                      fsg_hist_entry_rc(hist_entry));
            }
        }
    }
}


/*
 * Perform cross-word transitions; propagate each history entry created in this
 * frame to lextree roots attached to the target FSG state for that entry.
 */
static void
fsg_search_word_trans(fsg_search_t *fsgs)
{
    int32 bpidx, n_entries;
    fsg_hist_entry_t *hist_entry;
    fsg_link_t *l;
    int32 score, newscore, thresh, nf, d;
    fsg_pnode_t *root;
    int32 lc, rc;

    n_entries = fsg_history_n_entries(fsgs->history);

    thresh = fsgs->bestscore + fsgs->beam;
    nf = fsgs->frame + 1;

    for (bpidx = fsgs->bpidx_start; bpidx < n_entries; bpidx++) {
        hist_entry = fsg_history_entry_get(fsgs->history, bpidx);
        assert(hist_entry);
        score = fsg_hist_entry_score(hist_entry);
        assert(fsgs->frame == fsg_hist_entry_frame(hist_entry));

        l = fsg_hist_entry_fsglink(hist_entry);

        /* Destination state for hist_entry */
        d = l ? fsg_link_to_state(l) : fsg_model_start_state(fsgs->
                                                                fsg);

        lc = fsg_hist_entry_lc(hist_entry);

        /* Transition to all root nodes attached to state d */
        for (root = fsg_lextree_root(fsgs->lextree, d);
             root; root = root->sibling) {
            rc = root->ci_ext;

            if ((root->ctxt.bv[lc >> 5] & (1 << (lc & 0x001f))) &&
                (hist_entry->rc.bv[rc >> 5] & (1 << (rc & 0x001f)))) {
                /*
                 * Last CIphone of history entry is in left-context list supported by
                 * target root node, and
                 * first CIphone of target root node is in right context list supported
                 * by history entry;
                 * So the transition can go ahead (if new score is good enough).
                 */
                newscore = score + root->logs2prob;

                if ((newscore BETTER_THAN thresh)
                    && (newscore BETTER_THAN hmm_in_score(&root->hmm))) {
                    if (hmm_frame(&root->hmm) < nf) {
                        /* Newly activated node; add to active list */
                        fsgs->pnode_active_next =
                            glist_add_ptr(fsgs->pnode_active_next,
                                          (void *) root);
#if __FSG_DBG__
                        E_INFO
                            ("[%5d] WordTrans bpidx[%d] -> pnode[%08x] (activated)\n",
                             fsgs->frame, bpidx, (int32) root);
#endif
                    }
                    else {
#if __FSG_DBG__
                        E_INFO
                            ("[%5d] WordTrans bpidx[%d] -> pnode[%08x]\n",
                             fsgs->frame, bpidx, (int32) root);
#endif
                    }

                    hmm_enter(&root->hmm, newscore, bpidx, nf);
                }
            }
        }
    }
}


int
fsg_search_step(ps_search_t *search, int frame_idx)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;
    int16 const *senscr;
    acmod_t *acmod = search->acmod;
    gnode_t *gn;
    fsg_pnode_t *pnode;
    hmm_t *hmm;

    /* Activate our HMMs for the current frame if need be. */
    if (!acmod->compallsen)
        fsg_search_sen_active(fsgs);
    /* Compute GMM scores for the current frame. */
    senscr = acmod_score(acmod, &frame_idx);
    fsgs->n_sen_eval += acmod->n_senone_active;
    hmm_context_set_senscore(fsgs->hmmctx, senscr);

    /* Mark backpointer table for current frame. */
    fsgs->bpidx_start = fsg_history_n_entries(fsgs->history);

    /* Evaluate all active pnodes (HMMs) */
    fsg_search_hmm_eval(fsgs);

    /*
     * Prune and propagate the HMMs evaluated; create history entries for
     * word exits.  The words exits are tentative, and may be pruned; make
     * the survivors permanent via fsg_history_end_frame().
     */
    fsg_search_hmm_prune_prop(fsgs);
    fsg_history_end_frame(fsgs->history);

    /*
     * Propagate new history entries through any null transitions, creating
     * new history entries, and then make the survivors permanent.
     */
    fsg_search_null_prop(fsgs);
    fsg_history_end_frame(fsgs->history);

    /*
     * Perform cross-word transitions; propagate each history entry across its
     * terminating state to the root nodes of the lextree attached to the state.
     */
    fsg_search_word_trans(fsgs);

    /*
     * We've now come full circle, HMM and FSG states have been updated for
     * the next frame.
     * Update the active lists, deactivate any currently active HMMs that
     * did not survive into the next frame
     */
    for (gn = fsgs->pnode_active; gn; gn = gnode_next(gn)) {
        pnode = (fsg_pnode_t *) gnode_ptr(gn);
        hmm = fsg_pnode_hmmptr(pnode);

        if (hmm_frame(hmm) == fsgs->frame) {
            /* This HMM NOT activated for the next frame; reset it */
            fsg_psubtree_pnode_deactivate(pnode);
        }
        else {
            assert(hmm_frame(hmm) == (fsgs->frame + 1));
        }
    }

    /* Free the currently active list */
    glist_free(fsgs->pnode_active);

    /* Make the next-frame active list the current one */
    fsgs->pnode_active = fsgs->pnode_active_next;
    fsgs->pnode_active_next = NULL;

    /* End of this frame; ready for the next */
    ++fsgs->frame;

    return 1;
}


/*
 * Set all HMMs to inactive, clear active lists, initialize FSM start
 * state to be the only active node.
 * (Executed at the start of each utterance.)
 */
int
fsg_search_start(ps_search_t *search)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;
    int32 silcipid;
    fsg_pnode_ctxt_t ctxt;

    /* Reset dynamic adjustment factor for beams */
    fsgs->beam_factor = 1.0f;
    fsgs->beam = fsgs->beam_orig;
    fsgs->pbeam = fsgs->pbeam_orig;
    fsgs->wbeam = fsgs->wbeam_orig;

    silcipid = bin_mdef_ciphone_id(ps_search_acmod(fsgs)->mdef, "SIL");

    /* Initialize EVERYTHING to be inactive */
    assert(fsgs->pnode_active == NULL);
    assert(fsgs->pnode_active_next == NULL);

    fsg_history_reset(fsgs->history);
    fsg_history_utt_start(fsgs->history);
    fsgs->final = FALSE;

    /* Dummy context structure that allows all right contexts to use this entry */
    fsg_pnode_add_all_ctxt(&ctxt);

    /* Create dummy history entry leading to start state */
    fsgs->frame = -1;
    fsgs->bestscore = 0;
    fsg_history_entry_add(fsgs->history,
                          NULL, -1, 0, -1, silcipid, ctxt);
    fsgs->bpidx_start = 0;

    /* Propagate dummy history entry through NULL transitions from start state */
    fsg_search_null_prop(fsgs);

    /* Perform word transitions from this dummy history entry */
    fsg_search_word_trans(fsgs);

    /* Make the next-frame active list the current one */
    fsgs->pnode_active = fsgs->pnode_active_next;
    fsgs->pnode_active_next = NULL;

    ++fsgs->frame;

    fsgs->n_hmm_eval = 0;
    fsgs->n_sen_eval = 0;

    return 0;
}

/*
 * Cleanup at the end of each utterance.
 */
int
fsg_search_finish(ps_search_t *search)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;
    gnode_t *gn;
    fsg_pnode_t *pnode;
    int32 n_hist;

    /* Deactivate all nodes in the current and next-frame active lists */
    for (gn = fsgs->pnode_active; gn; gn = gnode_next(gn)) {
        pnode = (fsg_pnode_t *) gnode_ptr(gn);
        fsg_psubtree_pnode_deactivate(pnode);
    }
    for (gn = fsgs->pnode_active_next; gn; gn = gnode_next(gn)) {
        pnode = (fsg_pnode_t *) gnode_ptr(gn);
        fsg_psubtree_pnode_deactivate(pnode);
    }

    glist_free(fsgs->pnode_active);
    fsgs->pnode_active = NULL;
    glist_free(fsgs->pnode_active_next);
    fsgs->pnode_active_next = NULL;

    fsgs->final = TRUE;

    n_hist = fsg_history_n_entries(fsgs->history);
    E_INFO
        ("%d frames, %d HMMs (%d/fr), %d senones (%d/fr), %d history entries (%d/fr)\n\n",
         fsgs->frame, fsgs->n_hmm_eval,
         (fsgs->frame > 0) ? fsgs->n_hmm_eval / fsgs->frame : 0,
         fsgs->n_sen_eval,
         (fsgs->frame > 0) ? fsgs->n_sen_eval / fsgs->frame : 0,
         n_hist, (fsgs->frame > 0) ? n_hist / fsgs->frame : 0);

    return 0;
}

static int
fsg_search_find_exit(fsg_search_t *fsgs, int frame_idx, int final, int32 *out_score, int32* out_is_final)
{
    fsg_hist_entry_t *hist_entry;
    fsg_model_t *fsg;
    int bpidx, frm, last_frm, besthist;
    int32 bestscore;

    if (frame_idx == -1)
        frame_idx = fsgs->frame - 1;
    last_frm = frm = frame_idx;

    /* Scan backwards to find a word exit in frame_idx. */
    bpidx = fsg_history_n_entries(fsgs->history) - 1;
    while (bpidx > 0) {
        hist_entry = fsg_history_entry_get(fsgs->history, bpidx);
        if (fsg_hist_entry_frame(hist_entry) <= frame_idx) {
            frm = last_frm = fsg_hist_entry_frame(hist_entry);
            break;
        }
    }

    /* No hypothesis (yet). */
    if (bpidx <= 0) 
        return bpidx;

    /* Now find best word exit in this frame. */
    bestscore = INT_MIN;
    besthist = -1;
    fsg = fsgs->fsg;
    while (frm == last_frm) {
        fsg_link_t *fl;
        int32 score;

        fl = fsg_hist_entry_fsglink(hist_entry);
        score = fsg_hist_entry_score(hist_entry);
        
        if (fl == NULL)
	    break;

	/* Prefer final hypothesis */
	if (score == bestscore && fsg_link_to_state(fl) == fsg_model_final_state(fsg)) {
    	    besthist = bpidx;
	} else if (score BETTER_THAN bestscore) {
            /* Only enforce the final state constraint if this is a final hypothesis. */
            if ((!final)
                || fsg_link_to_state(fl) == fsg_model_final_state(fsg)) {
                bestscore = score;
                besthist = bpidx;
            }
        }
        
        --bpidx;
        if (bpidx < 0)
            break;
        hist_entry = fsg_history_entry_get(fsgs->history, bpidx);
        frm = fsg_hist_entry_frame(hist_entry);
    }

    /* Final state not reached. */
    if (besthist == -1) {
        E_ERROR("Final result does not match the grammar in frame %d\n", frame_idx);
        return -1;
    }

    /* This here's the one we want. */
    if (out_score)
        *out_score = bestscore;
    if (out_is_final) {
	fsg_link_t *fl;
	hist_entry = fsg_history_entry_get(fsgs->history, besthist);
	fl = fsg_hist_entry_fsglink(hist_entry);
	*out_is_final = (fsg_link_to_state(fl) == fsg_model_final_state(fsg));
    }
    return besthist;
}

/* FIXME: Mostly duplicated with ngram_search_bestpath(). */
static ps_latlink_t *
fsg_search_bestpath(ps_search_t *search, int32 *out_score, int backward)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;

    if (search->last_link == NULL) {
        search->last_link = ps_lattice_bestpath(search->dag, NULL,
                                                1.0, fsgs->ascale);
        if (search->last_link == NULL)
            return NULL;
        /* Also calculate betas so we can fill in the posterior
         * probability field in the segmentation. */
        if (search->post == 0)
            search->post = ps_lattice_posterior(search->dag, NULL, fsgs->ascale);
    }
    if (out_score)
        *out_score = search->last_link->path_scr + search->dag->final_node_ascr;
    return search->last_link;
}

char const *
fsg_search_hyp(ps_search_t *search, int32 *out_score, int32 *out_is_final)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;
    dict_t *dict = ps_search_dict(search);
    char *c;
    size_t len;
    int bp, bpidx;

    /* Get last backpointer table index. */
    bpidx = fsg_search_find_exit(fsgs, fsgs->frame, fsgs->final, out_score, out_is_final);
    /* No hypothesis (yet). */
    if (bpidx <= 0)
        return NULL;

    /* If bestpath is enabled and the utterance is complete, then run it. */
    if (fsgs->bestpath && fsgs->final) {
        ps_lattice_t *dag;
        ps_latlink_t *link;

        if ((dag = fsg_search_lattice(search)) == NULL) {
    	    E_WARN("Failed to obtain the lattice while bestpath enabled\n");
            return NULL;
        }
        if ((link = fsg_search_bestpath(search, out_score, FALSE)) == NULL) {
    	    E_WARN("Failed to find the bestpath in a lattice\n");
            return NULL;
        }
        return ps_lattice_hyp(dag, link);
    }

    bp = bpidx;
    len = 0;
    while (bp > 0) {
        fsg_hist_entry_t *hist_entry = fsg_history_entry_get(fsgs->history, bp);
        fsg_link_t *fl = fsg_hist_entry_fsglink(hist_entry);
        char const *baseword;
        int32 wid;

        bp = fsg_hist_entry_pred(hist_entry);
        wid = fsg_link_wid(fl);
        if (wid < 0 || fsg_model_is_filler(fsgs->fsg, wid))
            continue;
        baseword = dict_basestr(dict,
                                dict_wordid(dict,
                                            fsg_model_word_str(fsgs->fsg, wid)));
        len += strlen(baseword) + 1;
    }
    
    ckd_free(search->hyp_str);
    if (len == 0) {
	search->hyp_str = NULL;
	return search->hyp_str;
    }
    search->hyp_str = ckd_calloc(1, len);

    bp = bpidx;
    c = search->hyp_str + len - 1;
    while (bp > 0) {
        fsg_hist_entry_t *hist_entry = fsg_history_entry_get(fsgs->history, bp);
        fsg_link_t *fl = fsg_hist_entry_fsglink(hist_entry);
        char const *baseword;
        int32 wid;

        bp = fsg_hist_entry_pred(hist_entry);
        wid = fsg_link_wid(fl);
        if (wid < 0 || fsg_model_is_filler(fsgs->fsg, wid))
            continue;
        baseword = dict_basestr(dict,
                                dict_wordid(dict,
                                            fsg_model_word_str(fsgs->fsg, wid)));
        len = strlen(baseword);
        c -= len;
        memcpy(c, baseword, len);
        if (c > search->hyp_str) {
            --c;
            *c = ' ';
        }
    }

    return search->hyp_str;
}

static void
fsg_seg_bp2itor(ps_seg_t *seg, fsg_hist_entry_t *hist_entry)
{
    fsg_search_t *fsgs = (fsg_search_t *)seg->search;
    fsg_hist_entry_t *ph = NULL;
    int32 bp;

    if ((bp = fsg_hist_entry_pred(hist_entry)) >= 0)
        ph = fsg_history_entry_get(fsgs->history, bp);
    seg->word = fsg_model_word_str(fsgs->fsg, hist_entry->fsglink->wid);
    seg->ef = fsg_hist_entry_frame(hist_entry);
    seg->sf = ph ? fsg_hist_entry_frame(ph) + 1 : 0;
    /* This is kind of silly but it happens for null transitions. */
    if (seg->sf > seg->ef) seg->sf = seg->ef;
    seg->prob = 0; /* Bogus value... */
    /* "Language model" score = transition probability. */
    seg->lback = 1;
    seg->lscr = hist_entry->fsglink->logs2prob;
    if (ph) {
        /* FIXME: Not sure exactly how cross-word triphones are handled. */
        seg->ascr = hist_entry->score - ph->score - seg->lscr;
    }
    else
        seg->ascr = hist_entry->score - seg->lscr;
}


static void
fsg_seg_free(ps_seg_t *seg)
{
    fsg_seg_t *itor = (fsg_seg_t *)seg;
    ckd_free(itor->hist);
    ckd_free(itor);
}

static ps_seg_t *
fsg_seg_next(ps_seg_t *seg)
{
    fsg_seg_t *itor = (fsg_seg_t *)seg;

    if (++itor->cur == itor->n_hist) {
        fsg_seg_free(seg);
        return NULL;
    }

    fsg_seg_bp2itor(seg, itor->hist[itor->cur]);
    return seg;
}

static ps_segfuncs_t fsg_segfuncs = {
    /* seg_next */ fsg_seg_next,
    /* seg_free */ fsg_seg_free
};

static ps_seg_t *
fsg_search_seg_iter(ps_search_t *search, int32 *out_score)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;
    fsg_seg_t *itor;
    int bp, bpidx, cur;

    bpidx = fsg_search_find_exit(fsgs, fsgs->frame, fsgs->final, out_score, NULL);
    /* No hypothesis (yet). */
    if (bpidx <= 0)
        return NULL;

    /* If bestpath is enabled and the utterance is complete, then run it. */
    if (fsgs->bestpath && fsgs->final) {
        ps_lattice_t *dag;
        ps_latlink_t *link;

        if ((dag = fsg_search_lattice(search)) == NULL)
            return NULL;
        if ((link = fsg_search_bestpath(search, out_score, TRUE)) == NULL)
            return NULL;
        return ps_lattice_seg_iter(dag, link, 1.0);
    }

    /* Calling this an "iterator" is a bit of a misnomer since we have
     * to get the entire backtrace in order to produce it.  On the
     * other hand, all we actually need is the bptbl IDs, and we can
     * allocate a fixed-size array of them. */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &fsg_segfuncs;
    itor->base.search = search;
    itor->base.lwf = 1.0;
    itor->n_hist = 0;
    bp = bpidx;
    while (bp > 0) {
        fsg_hist_entry_t *hist_entry = fsg_history_entry_get(fsgs->history, bp);
        bp = fsg_hist_entry_pred(hist_entry);
        ++itor->n_hist;
    }
    if (itor->n_hist == 0) {
        ckd_free(itor);
        return NULL;
    }
    itor->hist = ckd_calloc(itor->n_hist, sizeof(*itor->hist));
    cur = itor->n_hist - 1;
    bp = bpidx;
    while (bp > 0) {
        fsg_hist_entry_t *hist_entry = fsg_history_entry_get(fsgs->history, bp);
        itor->hist[cur] = hist_entry;
        bp = fsg_hist_entry_pred(hist_entry);
        --cur;
    }

    /* Fill in relevant fields for first element. */
    fsg_seg_bp2itor((ps_seg_t *)itor, itor->hist[0]);
    
    return (ps_seg_t *)itor;
}

static int
fsg_search_prob(ps_search_t *search)
{
    fsg_search_t *fsgs = (fsg_search_t *)search;

    /* If bestpath is enabled and the utterance is complete, then run it. */
    if (fsgs->bestpath && fsgs->final) {
        ps_lattice_t *dag;
        ps_latlink_t *link;

        if ((dag = fsg_search_lattice(search)) == NULL)
            return 0;
        if ((link = fsg_search_bestpath(search, NULL, TRUE)) == NULL)
            return 0;
        return search->post;
    }
    else {
        /* FIXME: Give some kind of good estimate here, eventually. */
        return 0;
    }
}

static ps_latnode_t *
find_node(ps_lattice_t *dag, fsg_model_t *fsg, int sf, int32 wid, int32 node_id)
{
    ps_latnode_t *node;

    for (node = dag->nodes; node; node = node->next)
        if ((node->sf == sf) && (node->wid == wid) && (node->node_id == node_id))
            break;
    return node;
}

static ps_latnode_t *
new_node(ps_lattice_t *dag, fsg_model_t *fsg, int sf, int ef, int32 wid, int32 node_id, int32 ascr)
{
    ps_latnode_t *node;

    node = find_node(dag, fsg, sf, wid, node_id);

    if (node) {
        /* Update end frames. */
        if (node->lef == -1 || node->lef < ef)
            node->lef = ef;
        if (node->fef == -1 || node->fef > ef)
            node->fef = ef;
        /* Update best link score. */
        if (ascr BETTER_THAN node->info.best_exit)
            node->info.best_exit = ascr;
    }
    else {
        /* New node; link to head of list */
        node = listelem_malloc(dag->latnode_alloc);
        node->wid = wid;
        node->sf = sf;
        node->fef = node->lef = ef;
        node->reachable = FALSE;
        node->entries = NULL;
        node->exits = NULL;
        node->info.best_exit = ascr;
        node->node_id = node_id;

        node->next = dag->nodes;
        dag->nodes = node;
        ++dag->n_nodes;
    }

    return node;
}

static ps_latnode_t *
find_start_node_one(fsg_search_t *fsgs, ps_lattice_t *dag)
{
    ps_latnode_t *node;
    glist_t start = NULL;
    int nstart = 0;

    /* Look for all nodes starting in frame zero with some exits. */
    for (node = dag->nodes; node; node = node->next) {
        if (node->sf == 0 && node->exits) {
            E_INFO("Start node %s.%d:%d:%d\n",
                   fsg_model_word_str(fsgs->fsg, node->wid),
                   node->sf, node->fef, node->lef);
            start = glist_add_ptr(start, node);
            ++nstart;
        }
    }

    /* If there was more than one start node candidate, then we need
     * to create an artificial start node with epsilon transitions to
     * all of them. */
    if (nstart == 1) {
        node = gnode_ptr(start);
    }
    else {
        gnode_t *st;
        int wid;

        wid = fsg_model_word_add(fsgs->fsg, "<s>");
        if (fsgs->fsg->silwords)
            bitvec_set(fsgs->fsg->silwords, wid);
        node = new_node(dag, fsgs->fsg, 0, 0, wid, -1, 0);
        for (st = start; st; st = gnode_next(st))
            ps_lattice_link(dag, node, gnode_ptr(st), 0, 0);
    }
    glist_free(start);
    return node;
}

static ps_latnode_t *
find_end_node_one(fsg_search_t *fsgs, ps_lattice_t *dag)
{
    ps_latnode_t *node;
    glist_t end = NULL;
    int nend = 0;

    /* Look for all nodes ending in last frame with some entries. */
    for (node = dag->nodes; node; node = node->next) {
        if (node->lef == dag->n_frames - 1 && node->entries) {
            E_INFO("End node %s.%d:%d:%d (%d)\n",
                   fsg_model_word_str(fsgs->fsg, node->wid),
                   node->sf, node->fef, node->lef, node->info.best_exit);
            end = glist_add_ptr(end, node);
            ++nend;
        }
    }

    if (nend == 1) {
        node = gnode_ptr(end);
    }
    else if (nend == 0) {
        ps_latnode_t *last = NULL;
        int ef = 0;

        /* If there were no end node candidates, then just use the
         * node with the last exit frame. */
        for (node = dag->nodes; node; node = node->next) {
            if (node->lef > ef && node->entries) {
                last = node;
                ef = node->lef;
            }
        }
        node = last;
        if (node)
            E_INFO("End node %s.%d:%d:%d (%d)\n",
                   fsg_model_word_str(fsgs->fsg, node->wid),
                   node->sf, node->fef, node->lef, node->info.best_exit);
    }    
    else {
        /* If there was more than one end node candidate, then we need
         * to create an artificial end node with epsilon transitions
         * out of all of them. */
        gnode_t *st;
        int wid;
        wid = fsg_model_word_add(fsgs->fsg, "</s>");
        if (fsgs->fsg->silwords)
            bitvec_set(fsgs->fsg->silwords, wid);
        node = new_node(dag, fsgs->fsg, fsgs->frame, fsgs->frame, wid, -1, 0);
        /* Use the "best" (in reality it will be the only) exit link
         * score from this final node as the link score. */
        for (st = end; st; st = gnode_next(st)) {
            ps_latnode_t *src = gnode_ptr(st);
            ps_lattice_link(dag, src, node, src->info.best_exit, fsgs->frame);
        }
    }
    glist_free(end);
    return node;
}

static void
mark_reachable(ps_lattice_t *dag, ps_latnode_t *end)
{
    glist_t q;

    /* It doesn't matter which order we do this in. */
    end->reachable = TRUE;
    q = glist_add_ptr(NULL, end);
    while (q) {
        ps_latnode_t *node = gnode_ptr(q);
        latlink_list_t *x;

        /* Pop the front of the list. */
        q = gnode_free(q, NULL);
        /* Expand all its predecessors that haven't been seen yet. */
        for (x = node->entries; x; x = x->next) {
            ps_latnode_t *next = x->link->from;
            if (!next->reachable) {
                next->reachable = TRUE;
                q = glist_add_ptr(q, next);
            }
        }
    }
}

/**
 * Generate a lattice from FSG search results.
 *
 * One might think that this is simply a matter of adding acoustic
 * scores to the FSG's edges.  However, one would be wrong.  The
 * crucial difference here is that the word lattice is acyclic, and it
 * also contains timing information.
 */
static ps_lattice_t *
fsg_search_lattice(ps_search_t *search)
{
    fsg_search_t *fsgs;
    fsg_model_t *fsg;
    ps_latnode_t *node;
    ps_lattice_t *dag;
    int32 i, n;

    fsgs = (fsg_search_t *)search;

    /* Check to see if a lattice has previously been created over the
     * same number of frames, and reuse it if so. */
    if (search->dag && search->dag->n_frames == fsgs->frame)
        return search->dag;

    /* Nope, create a new one. */
    ps_lattice_free(search->dag);
    search->dag = NULL;
    dag = ps_lattice_init_search(search, fsgs->frame);
    fsg = fsgs->fsg;

    /*
     * Each history table entry represents a link in the word graph.
     * The set of nodes is determined by the number of unique
     * (word,start-frame) pairs in the history table.  So we will
     * first find all those nodes.
     */
    n = fsg_history_n_entries(fsgs->history);
    for (i = 0; i < n; ++i) {
        fsg_hist_entry_t *fh = fsg_history_entry_get(fsgs->history, i);
        int32 ascr;
        int sf;

        /* Skip null transitions. */
        if (fh->fsglink == NULL || fh->fsglink->wid == -1)
            continue;

        /* Find the start node of this link. */
        if (fh->pred) {
            fsg_hist_entry_t *pfh = fsg_history_entry_get(fsgs->history, fh->pred);
            /* FIXME: We include the transition score in the lattice
             * link score.  This is because of the practical
             * difficulty of obtaining it separately in bestpath or
             * forward-backward search, and because it is essentially
             * a unigram probability, so there is no need to treat it
             * separately from the acoustic score.  However, it's not
             * clear that this will actually yield correct results.*/
            ascr = fh->score - pfh->score;
            sf = pfh->frame + 1;
        }
        else {
            ascr = fh->score;
            sf = 0;
        }

        /*
         * Note that although scores are tied to links rather than
         * nodes, it's possible that there are no links out of the
         * destination node, and thus we need to preserve its score in
         * case it turns out to be utterance-final.
         */
        new_node(dag, fsg, sf, fh->frame, fh->fsglink->wid, fsg_link_to_state(fh->fsglink), ascr);
    }

    /*
     * Now, we will create links only to nodes that actually exist.
     */
    n = fsg_history_n_entries(fsgs->history);
    for (i = 0; i < n; ++i) {
        fsg_hist_entry_t *fh = fsg_history_entry_get(fsgs->history, i);
        fsg_arciter_t *itor;
        ps_latnode_t *src, *dest;
        int32 ascr;
        int sf;

        /* Skip null transitions. */
        if (fh->fsglink == NULL || fh->fsglink->wid == -1)
            continue;

        /* Find the start node of this link and calculate its link score. */
        if (fh->pred) {
            fsg_hist_entry_t *pfh = fsg_history_entry_get(fsgs->history, fh->pred);
            sf = pfh->frame + 1;
            ascr = fh->score - pfh->score;
        }
        else {
            ascr = fh->score;
            sf = 0;
        }
        src = find_node(dag, fsg, sf, fh->fsglink->wid, fsg_link_to_state(fh->fsglink));
        sf = fh->frame + 1;

        for (itor = fsg_model_arcs(fsg, fsg_link_to_state(fh->fsglink));
             itor; itor = fsg_arciter_next(itor)) {
            fsg_link_t *link = fsg_arciter_get(itor);
            
            /* FIXME: Need to figure out what to do about tag transitions. */
            if (link->wid >= 0) {
                /*
                 * For each non-epsilon link following this one, look for a
                 * matching node in the lattice and link to it.
                 */
                if ((dest = find_node(dag, fsg, sf, link->wid, fsg_link_to_state(link))) != NULL)
            	    ps_lattice_link(dag, src, dest, ascr, fh->frame);
            }
            else {
                /*
                 * Transitive closure on nulls has already been done, so we
                 * just need to look one link forward from them.
                 */
                fsg_arciter_t *itor2;
                
                /* Add all non-null links out of j. */
                for (itor2 = fsg_model_arcs(fsg, fsg_link_to_state(link));
                     itor2; itor2 = fsg_arciter_next(itor2)) {
                    fsg_link_t *link = fsg_arciter_get(itor2);

                    if (link->wid == -1)
                        continue;
                    
                    if ((dest = find_node(dag, fsg, sf, link->wid, fsg_link_to_state(link))) != NULL) {
                        ps_lattice_link(dag, src, dest, ascr, fh->frame);
                    }
                }
            }
        }
    }


    /* Figure out which nodes are the start and end nodes. */
    if ((dag->start = find_start_node_one(fsgs, dag)) == NULL) {
	E_WARN("Failed to find the start node\n");
        goto error_out;
    }
    if ((dag->end = find_end_node_one(fsgs, dag)) == NULL) {
	E_WARN("Failed to find the end node\n");
        goto error_out;
    }


    E_INFO("lattice start node %s.%d end node %s.%d\n",
           fsg_model_word_str(fsg, dag->start->wid), dag->start->sf,
           fsg_model_word_str(fsg, dag->end->wid), dag->end->sf);
    /* FIXME: Need to calculate final_node_ascr here. */

    /*
     * Convert word IDs from FSG to dictionary.
     */
    for (node = dag->nodes; node; node = node->next) {
        node->wid = dict_wordid(dag->search->dict,
                                fsg_model_word_str(fsg, node->wid));
        node->basewid = dict_basewid(dag->search->dict, node->wid);
    }

    /*
     * Now we are done, because the links in the graph are uniquely
     * defined by the history table.  However we should remove any
     * nodes which are not reachable from the end node of the FSG.
     * Everything is reachable from the start node by definition.
     */
    mark_reachable(dag, dag->end);

    ps_lattice_delete_unreachable(dag);
    {
        int32 silpen, fillpen;

        silpen = (int32)(logmath_log(fsg->lmath,
                                     cmd_ln_float32_r(ps_search_config(fsgs), "-silprob"))
                         * fsg->lw)
            >> SENSCR_SHIFT;
        fillpen = (int32)(logmath_log(fsg->lmath,
                                      cmd_ln_float32_r(ps_search_config(fsgs), "-fillprob"))
                          * fsg->lw)
            >> SENSCR_SHIFT;
        ps_lattice_bypass_fillers(dag, silpen, fillpen);
    }
    search->dag = dag;

    return dag;


error_out:
    ps_lattice_free(dag);
    return NULL;

}
/* -*- c-basic-offset:4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * fsg_search_internal.h -- Search structures for FSG decoding.
 */


#ifndef __S2_FSG_SEARCH_H__
#define __S2_FSG_SEARCH_H__


/* SphinxBase headers. */
#include <sphinxbase/glist.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/fsg_model.h>

/* Local headers. */
#include "pocketsphinx_internal.h"
#include "hmm.h"
#include "fsg_history.h"
#include "fsg_lextree.h"

/**
 * Segmentation "iterator" for FSG history.
 */
typedef struct fsg_seg_s {
    ps_seg_t base;  /**< Base structure. */
    fsg_hist_entry_t **hist;   /**< Sequence of history entries. */
    int16 n_hist;  /**< Number of history entries. */
    int16 cur;      /**< Current position in hist. */
} fsg_seg_t;

/**
 * Implementation of FSG search (and "FSG set") structure.
 */
typedef struct fsg_search_s {
    ps_search_t base;

    hmm_context_t *hmmctx; /**< HMM context. */

    hash_table_t *fsgs;		/**< Table of all FSGs loaded */
    fsg_model_t *fsg;		/**< Currently active FSG; NULL if none.  One
				   must be made active before starting FSG
				   decoding */
    jsgf_t *jsgf;               /**< Active JSGF grammar file. */
    struct fsg_lextree_s *lextree;/**< Lextree structure for the currently
				   active FSG */
    struct fsg_history_s *history;/**< For storing the Viterbi search history */
  
    glist_t pnode_active;	/**< Those active in this frame */
    glist_t pnode_active_next;	/**< Those activated for the next frame */
  
    int32 beam_orig;		/**< Global pruning threshold */
    int32 pbeam_orig;		/**< Pruning threshold for phone transition */
    int32 wbeam_orig;		/**< Pruning threshold for word exit */
    float32 beam_factor;	/**< Dynamic/adaptive factor (<=1) applied to above
                                     beams to determine actual effective beams.
                                     For implementing absolute pruning. */
    int32 beam, pbeam, wbeam;	/**< Effective beams after applying beam_factor */
    int32 lw, pip, wip;         /**< Language weights */
  
    frame_idx_t frame;		/**< Current frame. */
    uint8 final;		/**< Decoding is finished for this utterance. */
    uint8 bestpath;		/**< Whether to run bestpath search
                                   and confidence annotation at end. */
    float32 ascale;             /**< Acoustic score scale for posterior probabilities. */

    int32 bestscore;		/**< For beam pruning */
    int32 bpidx_start;		/**< First history entry index this frame */
  
    int32 ascr, lscr;		/**< Total acoustic and lm score for utt */
  
    int32 n_hmm_eval;		/**< Total HMMs evaluated this utt */
    int32 n_sen_eval;		/**< Total senones evaluated this utt */
} fsg_search_t;

/* Access macros */
#define fsg_search_frame(s)	((s)->frame)

/**
 * Create, initialize and return a search module.
 */
ps_search_t *fsg_search_init(cmd_ln_t *config,
                             acmod_t *acmod,
                             dict_t *dict,
                             dict2pid_t *d2p);

/**
 * Deallocate search structure.
 */
void fsg_search_free(ps_search_t *search);

/**
 * Update FSG search module for new or updated FSGs.
 */
int fsg_search_reinit(ps_search_t *fsgs, dict_t *dict, dict2pid_t *d2p);

/**
 * Prepare the FSG search structure for beginning decoding of the next
 * utterance.
 */
int fsg_search_start(ps_search_t *search);

/**
 * Step one frame forward through the Viterbi search.
 */
int fsg_search_step(ps_search_t *search, int frame_idx);

/**
 * Windup and clean the FSG search structure after utterance.
 */
int fsg_search_finish(ps_search_t *search);

/**
 * Get hypothesis string from the FSG search.
 */
char const *fsg_search_hyp(ps_search_t *search, int32 *out_score, int32 *out_is_final);

#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file hmm.h Implementation of HMM base structure.
 */

/* System headers. */
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "hmm.h"

hmm_context_t *
hmm_context_init(int32 n_emit_state,
		 uint8 ** const *tp,
		 int16 const *senscore,
		 uint16 * const *sseq)
{
    hmm_context_t *ctx;

    assert(n_emit_state > 0);
    if (n_emit_state > HMM_MAX_NSTATE) {
        E_ERROR("Number of emitting states must be <= %d\n", HMM_MAX_NSTATE);
        return NULL;
    }

    ctx = ckd_calloc(1, sizeof(*ctx));
    ctx->n_emit_state = n_emit_state;
    ctx->tp = tp;
    ctx->senscore = senscore;
    ctx->sseq = sseq;
    ctx->st_sen_scr = ckd_calloc(n_emit_state, sizeof(*ctx->st_sen_scr));

    return ctx;
}

void
hmm_context_free(hmm_context_t *ctx)
{
    if (ctx == NULL)
        return;
    ckd_free(ctx->st_sen_scr);
    ckd_free(ctx);
}

void
hmm_init(hmm_context_t *ctx, hmm_t *hmm, int mpx, int ssid, int tmatid)
{
    hmm->ctx = ctx;
    hmm->mpx = mpx;
    hmm->n_emit_state = ctx->n_emit_state;
    if (mpx) {
        int i;
        hmm->ssid = BAD_SSID;
        hmm->senid[0] = ssid;
        for (i = 1; i < hmm_n_emit_state(hmm); ++i) {
            hmm->senid[i] = BAD_SSID;
        }
    }
    else {
        hmm->ssid = ssid;
        memcpy(hmm->senid, ctx->sseq[ssid], hmm->n_emit_state * sizeof(*hmm->senid));
    }
    hmm->tmatid = tmatid;
    hmm_clear(hmm);
}

void
hmm_deinit(hmm_t *hmm)
{
}

void
hmm_dump(hmm_t * hmm,
         FILE * fp)
{
    int32 i;

    if (hmm_is_mpx(hmm)) {
        fprintf(fp, "MPX   ");
        for (i = 0; i < hmm_n_emit_state(hmm); i++)
            fprintf(fp, " %11d", hmm_senid(hmm, i));
        fprintf(fp, " ( ");
        for (i = 0; i < hmm_n_emit_state(hmm); i++)
            fprintf(fp, "%d ", hmm_ssid(hmm, i));
        fprintf(fp, ")\n");
    }
    else {
        fprintf(fp, "SSID  ");
        for (i = 0; i < hmm_n_emit_state(hmm); i++)
            fprintf(fp, " %11d", hmm_senid(hmm, i));
        fprintf(fp, " (%d)\n", hmm_ssid(hmm, 0));
    }

    if (hmm->ctx->senscore) {
        fprintf(fp, "SENSCR");
        for (i = 0; i < hmm_n_emit_state(hmm); i++)
            fprintf(fp, " %11d", hmm_senscr(hmm, i));
        fprintf(fp, "\n");
    }

    fprintf(fp, "SCORES %11d", hmm_in_score(hmm));
    for (i = 1; i < hmm_n_emit_state(hmm); i++)
        fprintf(fp, " %11d", hmm_score(hmm, i));
    fprintf(fp, " %11d", hmm_out_score(hmm));
    fprintf(fp, "\n");

    fprintf(fp, "HISTID %11d", hmm_in_history(hmm));
    for (i = 1; i < hmm_n_emit_state(hmm); i++)
        fprintf(fp, " %11d", hmm_history(hmm, i));
    fprintf(fp, " %11d", hmm_out_history(hmm));
    fprintf(fp, "\n");

    if (hmm_in_score(hmm) > 0)
        fprintf(fp,
                "ALERT!! The input score %d is large than 0. Probably wrap around.\n",
                hmm_in_score(hmm));
    if (hmm_out_score(hmm) > 0)
        fprintf(fp,
                "ALERT!! The output score %d is large than 0. Probably wrap around\n.",
                hmm_out_score(hmm));

    fflush(fp);
}


void
hmm_clear_scores(hmm_t * h)
{
    int32 i;

    hmm_in_score(h) = WORST_SCORE;
    for (i = 1; i < hmm_n_emit_state(h); i++)
        hmm_score(h, i) = WORST_SCORE;
    hmm_out_score(h) = WORST_SCORE;

    h->bestscore = WORST_SCORE;
}

void
hmm_clear(hmm_t * h)
{
    int32 i;

    hmm_in_score(h) = WORST_SCORE;
    hmm_in_history(h) = -1;
    for (i = 1; i < hmm_n_emit_state(h); i++) {
        hmm_score(h, i) = WORST_SCORE;
        hmm_history(h, i) = -1;
    }
    hmm_out_score(h) = WORST_SCORE;
    hmm_out_history(h) = -1;

    h->bestscore = WORST_SCORE;
    h->frame = -1;
}

void
hmm_enter(hmm_t *h, int32 score, int32 histid, int frame)
{
    hmm_in_score(h) = score;
    hmm_in_history(h) = histid;
    hmm_frame(h) = frame;
}

void
hmm_normalize(hmm_t *h, int32 bestscr)
{
    int32 i;

    for (i = 0; i < hmm_n_emit_state(h); i++) {
        if (hmm_score(h, i) BETTER_THAN WORST_SCORE)
            hmm_score(h, i) -= bestscr;
    }
    if (hmm_out_score(h) BETTER_THAN WORST_SCORE)
        hmm_out_score(h) -= bestscr;
}

#define hmm_tprob_5st(i, j) (-tp[(i)*6+(j)])
#define nonmpx_senscr(i) (-senscore[sseq[i]])

static int32
hmm_vit_eval_5st_lr(hmm_t * hmm)
{
    int16 const *senscore = hmm->ctx->senscore;
    uint8 const *tp = hmm->ctx->tp[hmm->tmatid][0];
    uint16 const *sseq = hmm->senid;
    int32 s5, s4, s3, s2, s1, s0, t2, t1, t0, bestScore;

    /* It was the best of scores, it was the worst of scores. */
    bestScore = WORST_SCORE;

    /* Cache problem here! */
    s4 = hmm_score(hmm, 4) + nonmpx_senscr(4);
    s3 = hmm_score(hmm, 3) + nonmpx_senscr(3);
    /* Transitions into non-emitting state 5 */
    if (s3 BETTER_THAN WORST_SCORE) {
        t1 = s4 + hmm_tprob_5st(4, 5);
        t2 = s3 + hmm_tprob_5st(3, 5);
        if (t1 BETTER_THAN t2) {
            s5 = t1;
            hmm_out_history(hmm)  = hmm_history(hmm, 4);
        } else {
            s5 = t2;
            hmm_out_history(hmm)  = hmm_history(hmm, 3);
        }
        if (s5 WORSE_THAN WORST_SCORE) s5 = WORST_SCORE;
        hmm_out_score(hmm) = s5;
        bestScore = s5;
    }

    s2 = hmm_score(hmm, 2) + nonmpx_senscr(2);
    /* All transitions into state 4 */
    if (s2 BETTER_THAN WORST_SCORE) {
        t0 = s4 + hmm_tprob_5st(4, 4);
        t1 = s3 + hmm_tprob_5st(3, 4);
        t2 = s2 + hmm_tprob_5st(2, 4);
        if (t0 BETTER_THAN t1) {
            if (t2 BETTER_THAN t0) {
                s4 = t2;
                hmm_history(hmm, 4)  = hmm_history(hmm, 2);
            } else
                s4 = t0;
        } else {
            if (t2 BETTER_THAN t1) {
                s4 = t2;
                hmm_history(hmm, 4)  = hmm_history(hmm, 2);
            } else {
                s4 = t1;
                hmm_history(hmm, 4)  = hmm_history(hmm, 3);
            }
        }
        if (s4 WORSE_THAN WORST_SCORE) s4 = WORST_SCORE;
        if (s4 BETTER_THAN bestScore) bestScore = s4;
        hmm_score(hmm, 4) = s4;
    }

    s1 = hmm_score(hmm, 1) + nonmpx_senscr(1);
    /* All transitions into state 3 */
    if (s1 BETTER_THAN WORST_SCORE) {
        t0 = s3 + hmm_tprob_5st(3, 3);
        t1 = s2 + hmm_tprob_5st(2, 3);
        t2 = s1 + hmm_tprob_5st(1, 3);
        if (t0 BETTER_THAN t1) {
            if (t2 BETTER_THAN t0) {
                s3 = t2;
                hmm_history(hmm, 3)  = hmm_history(hmm, 1);
            } else
                s3 = t0;
        } else {
            if (t2 BETTER_THAN t1) {
                s3 = t2;
                hmm_history(hmm, 3)  = hmm_history(hmm, 1);
            } else {
                s3 = t1;
                hmm_history(hmm, 3)  = hmm_history(hmm, 2);
            }
        }
        if (s3 WORSE_THAN WORST_SCORE) s3 = WORST_SCORE;
        if (s3 BETTER_THAN bestScore) bestScore = s3;
        hmm_score(hmm, 3) = s3;
    }

    s0 = hmm_in_score(hmm) + nonmpx_senscr(0);
    /* All transitions into state 2 (state 0 is always active) */
    t0 = s2 + hmm_tprob_5st(2, 2);
    t1 = s1 + hmm_tprob_5st(1, 2);
    t2 = s0 + hmm_tprob_5st(0, 2);
    if (t0 BETTER_THAN t1) {
        if (t2 BETTER_THAN t0) {
            s2 = t2;
            hmm_history(hmm, 2)  = hmm_in_history(hmm);
        } else
            s2 = t0;
    } else {
        if (t2 BETTER_THAN t1) {
            s2 = t2;
            hmm_history(hmm, 2)  = hmm_in_history(hmm);
        } else {
            s2 = t1;
            hmm_history(hmm, 2)  = hmm_history(hmm, 1);
        }
    }
    if (s2 WORSE_THAN WORST_SCORE) s2 = WORST_SCORE;
    if (s2 BETTER_THAN bestScore) bestScore = s2;
    hmm_score(hmm, 2) = s2;


    /* All transitions into state 1 */
    t0 = s1 + hmm_tprob_5st(1, 1);
    t1 = s0 + hmm_tprob_5st(0, 1);
    if (t0 BETTER_THAN t1) {
        s1 = t0;
    } else {
        s1 = t1;
        hmm_history(hmm, 1)  = hmm_in_history(hmm);
    }
    if (s1 WORSE_THAN WORST_SCORE) s1 = WORST_SCORE;
    if (s1 BETTER_THAN bestScore) bestScore = s1;
    hmm_score(hmm, 1) = s1;

    /* All transitions into state 0 */
    s0 = s0 + hmm_tprob_5st(0, 0);
    if (s0 WORSE_THAN WORST_SCORE) s0 = WORST_SCORE;
    if (s0 BETTER_THAN bestScore) bestScore = s0;
    hmm_in_score(hmm) = s0;

    hmm_bestscore(hmm) = bestScore;
    return bestScore;
}

#define mpx_senid(st) sseq[ssid[st]][st]
#define mpx_senscr(st) (-senscore[mpx_senid(st)])

static int32
hmm_vit_eval_5st_lr_mpx(hmm_t * hmm)
{
    uint8 const *tp = hmm->ctx->tp[hmm->tmatid][0];
    int16 const *senscore = hmm->ctx->senscore;
    uint16 * const *sseq = hmm->ctx->sseq;
    uint16 *ssid = hmm->senid;
    int32 bestScore;
    int32 s5, s4, s3, s2, s1, s0, t2, t1, t0;

    /* Don't propagate WORST_SCORE */
    if (ssid[4] == BAD_SSID)
        s4 = t1 = WORST_SCORE;
    else {
        s4 = hmm_score(hmm, 4) + mpx_senscr(4);
        t1 = s4 + hmm_tprob_5st(4, 5);
    }
    if (ssid[3] == BAD_SSID)
        s3 = t2 = WORST_SCORE;
    else {
        s3 = hmm_score(hmm, 3) + mpx_senscr(3);
        t2 = s3 + hmm_tprob_5st(3, 5);
    }
    if (t1 BETTER_THAN t2) {
        s5 = t1;
        hmm_out_history(hmm) = hmm_history(hmm, 4);
    }
    else {
        s5 = t2;
        hmm_out_history(hmm) = hmm_history(hmm, 3);
    }
    if (s5 WORSE_THAN WORST_SCORE) s5 = WORST_SCORE;
    hmm_out_score(hmm) = s5;
    bestScore = s5;

    /* Don't propagate WORST_SCORE */
    if (ssid[2] == BAD_SSID)
        s2 = t2 = WORST_SCORE;
    else {
        s2 = hmm_score(hmm, 2) + mpx_senscr(2);
        t2 = s2 + hmm_tprob_5st(2, 4);
    }

    t0 = t1 = WORST_SCORE;
    if (s4 != WORST_SCORE)
        t0 = s4 + hmm_tprob_5st(4, 4);
    if (s3 != WORST_SCORE)
        t1 = s3 + hmm_tprob_5st(3, 4);
    if (t0 BETTER_THAN t1) {
        if (t2 BETTER_THAN t0) {
            s4 = t2;
            hmm_history(hmm, 4) = hmm_history(hmm, 2);
            ssid[4] = ssid[2];
        }
        else
            s4 = t0;
    }
    else {
        if (t2 BETTER_THAN t1) {
            s4 = t2;
            hmm_history(hmm, 4) = hmm_history(hmm, 2);
            ssid[4] = ssid[2];
        }
        else {
            s4 = t1;
            hmm_history(hmm, 4) = hmm_history(hmm, 3);
            ssid[4] = ssid[3];
        }
    }
    if (s4 WORSE_THAN WORST_SCORE) s4 = WORST_SCORE;
    if (s4 BETTER_THAN bestScore)
        bestScore = s4;
    hmm_score(hmm, 4) = s4;

    /* Don't propagate WORST_SCORE */
    if (ssid[1] == BAD_SSID)
        s1 = t2 = WORST_SCORE;
    else {
        s1 = hmm_score(hmm, 1) + mpx_senscr(1);
        t2 = s1 + hmm_tprob_5st(1, 3);
    }
    t0 = t1 = WORST_SCORE;
    if (s3 != WORST_SCORE)
        t0 = s3 + hmm_tprob_5st(3, 3);
    if (s2 != WORST_SCORE)
        t1 = s2 + hmm_tprob_5st(2, 3);
    if (t0 BETTER_THAN t1) {
        if (t2 BETTER_THAN t0) {
            s3 = t2;
            hmm_history(hmm, 3) = hmm_history(hmm, 1);
            ssid[3] = ssid[1];
        }
        else
            s3 = t0;
    }
    else {
        if (t2 BETTER_THAN t1) {
            s3 = t2;
            hmm_history(hmm, 3) = hmm_history(hmm, 1);
            ssid[3] = ssid[1];
        }
        else {
            s3 = t1;
            hmm_history(hmm, 3) = hmm_history(hmm, 2);
            ssid[3] = ssid[2];
        }
    }
    if (s3 WORSE_THAN WORST_SCORE) s3 = WORST_SCORE;
    if (s3 BETTER_THAN bestScore) bestScore = s3;
    hmm_score(hmm, 3) = s3;

    /* State 0 is always active */
    s0 = hmm_in_score(hmm) + mpx_senscr(0);

    /* Don't propagate WORST_SCORE */
    t0 = t1 = WORST_SCORE;
    if (s2 != WORST_SCORE)
        t0 = s2 + hmm_tprob_5st(2, 2);
    if (s1 != WORST_SCORE)
        t1 = s1 + hmm_tprob_5st(1, 2);
    t2 = s0 + hmm_tprob_5st(0, 2);
    if (t0 BETTER_THAN t1) {
        if (t2 BETTER_THAN t0) {
            s2 = t2;
            hmm_history(hmm, 2) = hmm_in_history(hmm);
            ssid[2] = ssid[0];
        }
        else
            s2 = t0;
    }
    else {
        if (t2 BETTER_THAN t1) {
            s2 = t2;
            hmm_history(hmm, 2) = hmm_in_history(hmm);
            ssid[2] = ssid[0];
        }
        else {
            s2 = t1;
            hmm_history(hmm, 2) = hmm_history(hmm, 1);
            ssid[2] = ssid[1];
        }
    }
    if (s2 WORSE_THAN WORST_SCORE) s2 = WORST_SCORE;
    if (s2 BETTER_THAN bestScore) bestScore = s2;
    hmm_score(hmm, 2) = s2;

    /* Don't propagate WORST_SCORE */
    t0 = WORST_SCORE;
    if (s1 != WORST_SCORE)
        t0 = s1 + hmm_tprob_5st(1, 1);
    t1 = s0 + hmm_tprob_5st(0, 1);
    if (t0 BETTER_THAN t1) {
        s1 = t0;
    }
    else {
        s1 = t1;
        hmm_history(hmm, 1) = hmm_in_history(hmm);
        ssid[1] = ssid[0];
    }
    if (s1 WORSE_THAN WORST_SCORE) s1 = WORST_SCORE;
    if (s1 BETTER_THAN bestScore) bestScore = s1;
    hmm_score(hmm, 1) = s1;

    s0 += hmm_tprob_5st(0, 0);
    if (s0 WORSE_THAN WORST_SCORE) s0 = WORST_SCORE;
    if (s0 BETTER_THAN bestScore) bestScore = s0;
    hmm_in_score(hmm) = s0;

    hmm_bestscore(hmm) = bestScore;
    return bestScore;
}

#define hmm_tprob_3st(i, j) (-tp[(i)*4+(j)])

static int32
hmm_vit_eval_3st_lr(hmm_t * hmm)
{
    int16 const *senscore = hmm->ctx->senscore;
    uint8 const *tp = hmm->ctx->tp[hmm->tmatid][0];
    uint16 const *sseq = hmm->senid;
    int32 s3, s2, s1, s0, t2, t1, t0, bestScore;

    s2 = hmm_score(hmm, 2) + nonmpx_senscr(2);
    s1 = hmm_score(hmm, 1) + nonmpx_senscr(1);
    s0 = hmm_in_score(hmm) + nonmpx_senscr(0);

    /* It was the best of scores, it was the worst of scores. */
    bestScore = WORST_SCORE;
    t2 = INT_MIN; /* Not used unless skipstate is true */

    /* Transitions into non-emitting state 3 */
    if (s1 BETTER_THAN WORST_SCORE) {
        t1 = s2 + hmm_tprob_3st(2, 3);
        if (hmm_tprob_3st(1,3) BETTER_THAN TMAT_WORST_SCORE)
            t2 = s1 + hmm_tprob_3st(1, 3);
        if (t1 BETTER_THAN t2) {
            s3 = t1;
            hmm_out_history(hmm)  = hmm_history(hmm, 2);
        } else {
            s3 = t2;
            hmm_out_history(hmm)  = hmm_history(hmm, 1);
        }
        if (s3 WORSE_THAN WORST_SCORE) s3 = WORST_SCORE;
        hmm_out_score(hmm) = s3;
        bestScore = s3;
    }

    /* All transitions into state 2 (state 0 is always active) */
    t0 = s2 + hmm_tprob_3st(2, 2);
    t1 = s1 + hmm_tprob_3st(1, 2);
    if (hmm_tprob_3st(0, 2) BETTER_THAN TMAT_WORST_SCORE)
        t2 = s0 + hmm_tprob_3st(0, 2);
    if (t0 BETTER_THAN t1) {
        if (t2 BETTER_THAN t0) {
            s2 = t2;
            hmm_history(hmm, 2)  = hmm_in_history(hmm);
        } else
            s2 = t0;
    } else {
        if (t2 BETTER_THAN t1) {
            s2 = t2;
            hmm_history(hmm, 2)  = hmm_in_history(hmm);
        } else {
            s2 = t1;
            hmm_history(hmm, 2)  = hmm_history(hmm, 1);
        }
    }
    if (s2 WORSE_THAN WORST_SCORE) s2 = WORST_SCORE;
    if (s2 BETTER_THAN bestScore) bestScore = s2;
    hmm_score(hmm, 2) = s2;

    /* All transitions into state 1 */
    t0 = s1 + hmm_tprob_3st(1, 1);
    t1 = s0 + hmm_tprob_3st(0, 1);
    if (t0 BETTER_THAN t1) {
        s1 = t0;
    } else {
        s1 = t1;
        hmm_history(hmm, 1)  = hmm_in_history(hmm);
    }
    if (s1 WORSE_THAN WORST_SCORE) s1 = WORST_SCORE;
    if (s1 BETTER_THAN bestScore) bestScore = s1;
    hmm_score(hmm, 1) = s1;

    /* All transitions into state 0 */
    s0 = s0 + hmm_tprob_3st(0, 0);
    if (s0 WORSE_THAN WORST_SCORE) s0 = WORST_SCORE;
    if (s0 BETTER_THAN bestScore) bestScore = s0;
    hmm_in_score(hmm) = s0;

    hmm_bestscore(hmm) = bestScore;
    return bestScore;
}

static int32
hmm_vit_eval_3st_lr_mpx(hmm_t * hmm)
{
    uint8 const *tp = hmm->ctx->tp[hmm->tmatid][0];
    int16 const *senscore = hmm->ctx->senscore;
    uint16 * const *sseq = hmm->ctx->sseq;
    uint16 *ssid = hmm->senid;
    int32 bestScore;
    int32 s3, s2, s1, s0, t2, t1, t0;

    /* Don't propagate WORST_SCORE */
    t2 = INT_MIN; /* Not used unless skipstate is true */
    if (ssid[2] == BAD_SSID)
        s2 = t1 = WORST_SCORE;
    else {
        s2 = hmm_score(hmm, 2) + mpx_senscr(2);
        t1 = s2 + hmm_tprob_3st(2, 3);
    }
    if (ssid[1] == BAD_SSID)
        s1 = t2 = WORST_SCORE;
    else {
        s1 = hmm_score(hmm, 1) + mpx_senscr(1);
        if (hmm_tprob_3st(1,3) BETTER_THAN TMAT_WORST_SCORE)
            t2 = s1 + hmm_tprob_3st(1, 3);
    }
    if (t1 BETTER_THAN t2) {
        s3 = t1;
        hmm_out_history(hmm) = hmm_history(hmm, 2);
    }
    else {
        s3 = t2;
        hmm_out_history(hmm) = hmm_history(hmm, 1);
    }
    if (s3 WORSE_THAN WORST_SCORE) s3 = WORST_SCORE;
    hmm_out_score(hmm) = s3;
    bestScore = s3;

    /* State 0 is always active */
    s0 = hmm_in_score(hmm) + mpx_senscr(0);

    /* Don't propagate WORST_SCORE */
    t0 = t1 = WORST_SCORE;
    if (s2 != WORST_SCORE)
        t0 = s2 + hmm_tprob_3st(2, 2);
    if (s1 != WORST_SCORE)
        t1 = s1 + hmm_tprob_3st(1, 2);
    if (hmm_tprob_3st(0,2) BETTER_THAN TMAT_WORST_SCORE)
        t2 = s0 + hmm_tprob_3st(0, 2);
    if (t0 BETTER_THAN t1) {
        if (t2 BETTER_THAN t0) {
            s2 = t2;
            hmm_history(hmm, 2) = hmm_in_history(hmm);
            ssid[2] = ssid[0];
        }
        else
            s2 = t0;
    }
    else {
        if (t2 BETTER_THAN t1) {
            s2 = t2;
            hmm_history(hmm, 2) = hmm_in_history(hmm);
            ssid[2] = ssid[0];
        }
        else {
            s2 = t1;
            hmm_history(hmm, 2) = hmm_history(hmm, 1);
            ssid[2] = ssid[1];
        }
    }
    if (s2 WORSE_THAN WORST_SCORE) s2 = WORST_SCORE;
    if (s2 BETTER_THAN bestScore) bestScore = s2;
    hmm_score(hmm, 2) = s2;

    /* Don't propagate WORST_SCORE */
    t0 = WORST_SCORE;
    if (s1 != WORST_SCORE)
        t0 = s1 + hmm_tprob_3st(1, 1);
    t1 = s0 + hmm_tprob_3st(0, 1);
    if (t0 BETTER_THAN t1) {
        s1 = t0;
    }
    else {
        s1 = t1;
        hmm_history(hmm, 1) = hmm_in_history(hmm);
        ssid[1] = ssid[0];
    }
    if (s1 WORSE_THAN WORST_SCORE) s1 = WORST_SCORE;
    if (s1 BETTER_THAN bestScore) bestScore = s1;
    hmm_score(hmm, 1) = s1;

    /* State 0 is always active */
    s0 += hmm_tprob_3st(0, 0);
    if (s0 WORSE_THAN WORST_SCORE) s0 = WORST_SCORE;
    if (s0 BETTER_THAN bestScore) bestScore = s0;
    hmm_in_score(hmm) = s0;

    hmm_bestscore(hmm) = bestScore;
    return bestScore;
}

static int32
hmm_vit_eval_anytopo(hmm_t * hmm)
{
    hmm_context_t *ctx = hmm->ctx;
    int32 to, from, bestfrom;
    int32 newscr, scr, bestscr;
    int final_state;

    /* Compute previous state-score + observation output prob for each emitting state */
    ctx->st_sen_scr[0] = hmm_in_score(hmm) + hmm_senscr(hmm, 0);
    for (from = 1; from < hmm_n_emit_state(hmm); ++from) {
        if ((ctx->st_sen_scr[from] =
             hmm_score(hmm, from) + hmm_senscr(hmm, from)) WORSE_THAN WORST_SCORE)
            ctx->st_sen_scr[from] = WORST_SCORE;
    }

    /* FIXME/TODO: Use the BLAS for all this. */
    /* Evaluate final-state first, which does not have a self-transition */
    final_state = hmm_n_emit_state(hmm);
    to = final_state;
    scr = WORST_SCORE;
    bestfrom = -1;
    for (from = to - 1; from >= 0; --from) {
        if ((hmm_tprob(hmm, from, to) BETTER_THAN TMAT_WORST_SCORE) &&
            ((newscr = ctx->st_sen_scr[from]
              + hmm_tprob(hmm, from, to)) BETTER_THAN scr)) {
            scr = newscr;
            bestfrom = from;
        }
    }
    hmm_out_score(hmm) = scr;
    if (bestfrom >= 0)
        hmm_out_history(hmm) = hmm_history(hmm, bestfrom);
    bestscr = scr;

    /* Evaluate all other states, which might have self-transitions */
    for (to = final_state - 1; to >= 0; --to) {
        /* Score from self-transition, if any */
        scr =
            (hmm_tprob(hmm, to, to) BETTER_THAN TMAT_WORST_SCORE)
            ? ctx->st_sen_scr[to] + hmm_tprob(hmm, to, to)
            : WORST_SCORE;

        /* Scores from transitions from other states */
        bestfrom = -1;
        for (from = to - 1; from >= 0; --from) {
            if ((hmm_tprob(hmm, from, to) BETTER_THAN TMAT_WORST_SCORE) &&
                ((newscr = ctx->st_sen_scr[from]
                  + hmm_tprob(hmm, from, to)) BETTER_THAN scr)) {
                scr = newscr;
                bestfrom = from;
            }
        }

        /* Update new result for state to */
        if (to == 0) {
            hmm_in_score(hmm) = scr;
            if (bestfrom >= 0)
                hmm_in_history(hmm) = hmm_history(hmm, bestfrom);
        }
        else {
            hmm_score(hmm, to) = scr;
            if (bestfrom >= 0)
                hmm_history(hmm, to) = hmm_history(hmm, bestfrom);
        }
        /* Propagate ssid for multiplex HMMs */
        if (bestfrom >= 0 && hmm_is_mpx(hmm))
            hmm->senid[to] = hmm->senid[bestfrom];

        if (bestscr WORSE_THAN scr)
            bestscr = scr;
    }

    hmm_bestscore(hmm) = bestscr;
    return bestscr;
}

int32
hmm_vit_eval(hmm_t * hmm)
{
    if (hmm_is_mpx(hmm)) {
        if (hmm_n_emit_state(hmm) == 5)
            return hmm_vit_eval_5st_lr_mpx(hmm);
        else if (hmm_n_emit_state(hmm) == 3)
            return hmm_vit_eval_3st_lr_mpx(hmm);
        else
            return hmm_vit_eval_anytopo(hmm);
    }
    else {
        if (hmm_n_emit_state(hmm) == 5)
            return hmm_vit_eval_5st_lr(hmm);
        else if (hmm_n_emit_state(hmm) == 3)
            return hmm_vit_eval_3st_lr(hmm);
        else
            return hmm_vit_eval_anytopo(hmm);
    }
}

int32
hmm_dump_vit_eval(hmm_t * hmm, FILE * fp)
{
    int32 bs = 0;

    if (fp) {
        fprintf(fp, "BEFORE:\n");
        hmm_dump(hmm, fp);
    }
    bs = hmm_vit_eval(hmm);
    if (fp) {
        fprintf(fp, "AFTER:\n");
        hmm_dump(hmm, fp);
    }

    return bs;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file hmm.h Hidden Markov Model base structures.
 */

#ifndef __HMM_H__
#define __HMM_H__

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/fixpoint.h>
#include <sphinxbase/listelem_alloc.h>

/* PocketSphinx headers. */
#include "bin_mdef.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

/** 
 * Type for frame index values. Used in HMM indexes and 
 * backpointers and affects memory required. Make it int32 to be able to 
 * process longer utterances. Due to limitations of FSG search implementation
 * this value needs to be signed.
 */
typedef int16 frame_idx_t;

/**
 * Maximum number of frames in index, should be in sync with above.
 */
#define MAX_N_FRAMES MAX_INT16


/** Shift count for senone scores. */
#define SENSCR_SHIFT 10

/**
 * Large "bad" score.
 *
 * This number must be "bad" enough so that 4 times WORST_SCORE will
 * not overflow. The reason for this is that the search doesn't check
 * the scores in a model before evaluating the model and it may
 * require as many was 4 plies before the new 'good' score can wipe
 * out the initial WORST_SCORE initialization.
 */
#define WORST_SCORE		((int)0xE0000000)

/**
 * Watch out, though!  Transition matrix entries that are supposed to
 * be "zero" don't actually get that small due to quantization.
 */
#define TMAT_WORST_SCORE	(-255)

/**
 * Is one score better than another?
 */
#define BETTER_THAN >

/**
 * Is one score worse than another?
 */
#define WORSE_THAN <

/** \file hmm.h
 * \brief HMM data structure and operation
 *
 * For efficiency, this version is hardwired for two possible HMM
 * topologies, but will fall back to others:
 * 
 * 5-state left-to-right HMMs: (0 is the *emitting* entry state and E
 * is a non-emitting exit state; the x's indicate allowed transitions
 * between source and destination states):
 * 
 * <pre>
 *               0   1   2   3   4   E (destination-states)
 *           0   x   x   x
 *           1       x   x   x
 *           2           x   x   x
 *           3               x   x   x
 *           4                   x   x
 *    (source-states)
 * </pre>
 *
 * 5-state topologies that contain a subset of the above transitions should work as well.
 * 
 * 3-state left-to-right HMMs (similar notation as the 5-state topology above):
 * 
 * <pre>
 *               0   1   2   E (destination-states)
 *           0   x   x   x
 *           1       x   x   x
 *           2           x   x 
 *    (source-states)
 * </pre>
 *
 * 3-state topologies that contain a subset of the above transitions should work as well. 
 */

/**
 * @struct hmm_context_t
 * @brief Shared information between a set of HMMs.
 *
 * We assume that the initial state is emitting and that the
 * transition matrix is n_emit_state x (n_emit_state+1), where the
 * extra destination dimension correponds to the non-emitting final or
 * exit state.
 */
typedef struct hmm_context_s {
    int32 n_emit_state;     /**< Number of emitting states in this set of HMMs. */
    uint8 ** const *tp;	    /**< State transition scores tp[id][from][to] (logs3 values). */
    int16 const *senscore;  /**< State emission scores senscore[senid]
                               (negated scaled logs3 values). */
    uint16 * const *sseq;   /**< Senone sequence mapping. */
    int32 *st_sen_scr;      /**< Temporary array of senone scores (for some topologies). */
    listelem_alloc_t *mpx_ssid_alloc; /**< Allocator for senone sequence ID arrays. */
    void *udata;            /**< Whatever you feel like, gosh. */
} hmm_context_t;

/**
 * Hard-coded limit on the number of emitting states.
 */
#define HMM_MAX_NSTATE 5

/**
 * @struct hmm_t
 * @brief An individual HMM among the HMM search space.
 *
 * An individual HMM among the HMM search space.  An HMM with N
 * emitting states consists of N+1 internal states including the
 * non-emitting exit (out) state.
 */
typedef struct hmm_s {
    hmm_context_t *ctx;            /**< Shared context data for this HMM. */
    int32 score[HMM_MAX_NSTATE];   /**< State scores for emitting states. */
    int32 history[HMM_MAX_NSTATE]; /**< History indices for emitting states. */
    int32 out_score;               /**< Score for non-emitting exit state. */
    int32 out_history;             /**< History index for non-emitting exit state. */
    uint16 ssid;                   /**< Senone sequence ID (for non-MPX) */
    uint16 senid[HMM_MAX_NSTATE];  /**< Senone IDs (non-MPX) or sequence IDs (MPX) */
    int32 bestscore;	/**< Best [emitting] state score in current frame (for pruning). */
    int16 tmatid;       /**< Transition matrix ID (see hmm_context_t). */
    frame_idx_t frame;  /**< Frame in which this HMM was last active; <0 if inactive */
    uint8 mpx;          /**< Is this HMM multiplex? (hoisted for speed) */
    uint8 n_emit_state; /**< Number of emitting states (hoisted for speed) */
} hmm_t;

/** Access macros. */
#define hmm_context(h) (h)->ctx
#define hmm_is_mpx(h) (h)->mpx

#define hmm_in_score(h) (h)->score[0]
#define hmm_score(h,st) (h)->score[st]
#define hmm_out_score(h) (h)->out_score

#define hmm_in_history(h) (h)->history[0]
#define hmm_history(h,st) (h)->history[st]
#define hmm_out_history(h) (h)->out_history

#define hmm_bestscore(h) (h)->bestscore
#define hmm_frame(h) (h)->frame
#define hmm_mpx_ssid(h,st) (h)->senid[st]
#define hmm_nonmpx_ssid(h) (h)->ssid
#define hmm_ssid(h,st) (hmm_is_mpx(h)                                   \
                        ? hmm_mpx_ssid(h,st) : hmm_nonmpx_ssid(h))
#define hmm_mpx_senid(h,st) (hmm_mpx_ssid(h,st) == BAD_SENID \
                             ? BAD_SENID : (h)->ctx->sseq[hmm_mpx_ssid(h,st)][st])
#define hmm_nonmpx_senid(h,st) ((h)->senid[st])
#define hmm_senid(h,st) (hmm_is_mpx(h)                                  \
                         ? hmm_mpx_senid(h,st) : hmm_nonmpx_senid(h,st))
#define hmm_senscr(h,st) (hmm_senid(h,st) == BAD_SENID                  \
                          ? WORST_SCORE                                 \
                          : -(h)->ctx->senscore[hmm_senid(h,st)])
#define hmm_tmatid(h) (h)->tmatid
#define hmm_tprob(h,i,j) (-(h)->ctx->tp[hmm_tmatid(h)][i][j])
#define hmm_n_emit_state(h) ((h)->n_emit_state)
#define hmm_n_state(h) ((h)->n_emit_state + 1)

/**
 * Create an HMM context.
 **/
hmm_context_t *hmm_context_init(int32 n_emit_state,
                                uint8 ** const *tp,
                                int16 const *senscore,
                                uint16 * const *sseq);

/**
 * Change the senone score array for a context.
 **/
#define hmm_context_set_senscore(ctx, senscr) ((ctx)->senscore = (senscr))

/**
 * Free an HMM context.
 *
 * @note The transition matrices, senone scores, and senone sequence
 * mapping are all assumed to be allocated externally, and will NOT be
 * freed by this function.
 **/
void hmm_context_free(hmm_context_t *ctx);

/**
 * Populate a previously-allocated HMM structure, allocating internal data.
 **/
void hmm_init(hmm_context_t *ctx, hmm_t *hmm, int mpx, int ssid, int tmatid);

/**
 * Free an HMM structure, releasing internal data (but not the HMM structure itself).
 */
void hmm_deinit(hmm_t *hmm);

/**
 * Reset the states of the HMM to the invalid condition.

 * i.e., scores to WORST_SCORE and hist to undefined.
 */
void hmm_clear(hmm_t *h);

/**
 * Reset the scores of the HMM.
 */
void hmm_clear_scores(hmm_t *h);

/**
 * Renormalize the scores in this HMM based on the given best score.
 */
void hmm_normalize(hmm_t *h, int32 bestscr);

/**
 * Enter an HMM with the given path score and history ID.
 **/
void hmm_enter(hmm_t *h, int32 score,
               int32 histid, int frame);

/**
 * Viterbi evaluation of given HMM.
 *
 * @note If this module were being used for tracking state
 * segmentations, the dummy, non-emitting exit state would have to be
 * updated separately.  In the Viterbi DP diagram, transitions to the
 * exit state occur from the current time; they are vertical
 * transitions.  Hence they should be made only after the history has
 * been logged for the emitting states.  But we're not bothered with
 * state segmentations, for now.  So, we update the exit state as
 * well.
*/
int32 hmm_vit_eval(hmm_t *hmm);
  

/**
 * Like hmm_vit_eval, but dump HMM state and relevant senscr to fp first, for debugging;.
 */
int32 hmm_dump_vit_eval(hmm_t *hmm,  /**< In/Out: HMM being updated */
                        FILE *fp /**< An output file pointer */
    );

/** 
 * For debugging, dump the whole HMM out.
 */

void hmm_dump(hmm_t *h,  /**< In/Out: HMM being updated */
              FILE *fp /**< An output file pointer */
    );


#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif

#endif /* __HMM_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * mdef.c -- HMM model definition: base (CI) phones and triphones
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1999 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * 
 * 
 * 22-Nov-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Imported from s3.2, for supporting s3 format continuous
 * 		acoustic models.
 * 
 * 14-Oct-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon
 * 		Added mdef_sseq2sen_active().
 * 
 * 06-May-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon
 * 		In mdef_phone_id(), added backing off to silence phone context from filler
 * 		context if original triphone not found.
 * 
 * 30-Apr-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon
 * 		Added senone-sequence id (ssid) to phone_t and appropriate functions to
 * 		maintain it.  Instead, moved state sequence info to mdef_t.
 * 
 * 13-Jul-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added mdef_phone_str().
 * 
 * 01-Jan-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Allowed mdef_phone_id_nearest to return base phone id if either
 * 		left or right context (or both) is undefined.
 * 
 * 01-Jan-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Created.
 */


/*
 * Major assumptions:
 *   All phones have same #states, same topology.
 *   Every phone has exactly one non-emitting, final state--the last one.
 *   CI phones must appear first in model definition file.
 */

/* System headers. */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "mdef.h"


#define MODEL_DEF_VERSION	"0.3"

static void
ciphone_add(mdef_t * m, char *ci, int p)
{
    assert(p < m->n_ciphone);

    m->ciphone[p].name = (char *) ckd_salloc(ci);       /* freed in mdef_free */
    if (hash_table_enter(m->ciphone_ht, m->ciphone[p].name,
                         (void *)(long)p) != (void *)(long)p)
        E_FATAL("hash_table_enter(%s) failed; duplicate CIphone?\n",
                m->ciphone[p].name);
}


static ph_lc_t *
find_ph_lc(ph_lc_t * lclist, int lc)
{
    ph_lc_t *lcptr;

    for (lcptr = lclist; lcptr && (lcptr->lc != lc); lcptr = lcptr->next);
    return lcptr;
}


static ph_rc_t *
find_ph_rc(ph_rc_t * rclist, int rc)
{
    ph_rc_t *rcptr;

    for (rcptr = rclist; rcptr && (rcptr->rc != rc); rcptr = rcptr->next);
    return rcptr;
}


static void
triphone_add(mdef_t * m,
             int ci, int lc, int rc, word_posn_t wpos,
             int p)
{
    ph_lc_t *lcptr;
    ph_rc_t *rcptr;

    assert(p < m->n_phone);

    /* Fill in phone[p] information (state and tmat mappings added later) */
    m->phone[p].ci = ci;
    m->phone[p].lc = lc;
    m->phone[p].rc = rc;
    m->phone[p].wpos = wpos;

    /* Create <ci,lc,rc,wpos> -> p mapping if not a CI phone */
    if (p >= m->n_ciphone) {
        if ((lcptr = find_ph_lc(m->wpos_ci_lclist[wpos][(int) ci], lc))
            == NULL) {
            lcptr = (ph_lc_t *) ckd_calloc(1, sizeof(ph_lc_t)); /* freed at mdef_free, I believe */
            lcptr->lc = lc;
            lcptr->next = m->wpos_ci_lclist[wpos][(int) ci];
            m->wpos_ci_lclist[wpos][(int) ci] = lcptr;  /* This is what needs to be freed */
        }
        if ((rcptr = find_ph_rc(lcptr->rclist, rc)) != NULL) {
            __BIGSTACKVARIABLE__ char buf[4096];

            mdef_phone_str(m, rcptr->pid, buf);
            E_FATAL("Duplicate triphone: %s\n", buf);
        }

        rcptr = (ph_rc_t *) ckd_calloc(1, sizeof(ph_rc_t));     /* freed in mdef_free, I believe */
        rcptr->rc = rc;
        rcptr->pid = p;
        rcptr->next = lcptr->rclist;
        lcptr->rclist = rcptr;
    }
}


int
mdef_ciphone_id(mdef_t * m, char *ci)
{
    int32 id;
    if (hash_table_lookup_int32(m->ciphone_ht, ci, &id) < 0)
        return -1;
    return id;
}


const char *
mdef_ciphone_str(mdef_t * m, int id)
{
    assert(m);
    assert((id >= 0) && (id < m->n_ciphone));

    return (m->ciphone[id].name);
}


int
mdef_phone_str(mdef_t * m, int pid, char *buf)
{
    char *wpos_name;

    assert(m);
    assert((pid >= 0) && (pid < m->n_phone));
    wpos_name = WPOS_NAME;

    buf[0] = '\0';
    if (pid < m->n_ciphone)
        sprintf(buf, "%s", mdef_ciphone_str(m, pid));
    else {
        sprintf(buf, "%s %s %s %c",
                mdef_ciphone_str(m, m->phone[pid].ci),
                mdef_ciphone_str(m, m->phone[pid].lc),
                mdef_ciphone_str(m, m->phone[pid].rc),
                wpos_name[m->phone[pid].wpos]);
    }
    return 0;
}


int
mdef_phone_id(mdef_t * m,
              int ci, int lc, int rc, word_posn_t wpos)
{
    ph_lc_t *lcptr;
    ph_rc_t *rcptr;
    int newl, newr;

    assert(m);
    assert((ci >= 0) && (ci < m->n_ciphone));
    assert((lc >= 0) && (lc < m->n_ciphone));
    assert((rc >= 0) && (rc < m->n_ciphone));
    assert((wpos >= 0) && (wpos < N_WORD_POSN));

    if (((lcptr =
          find_ph_lc(m->wpos_ci_lclist[wpos][(int) ci], lc)) == NULL)
        || ((rcptr = find_ph_rc(lcptr->rclist, rc)) == NULL)) {
        /* Not found; backoff to silence context if non-silence filler context */
        if (m->sil < 0)
            return -1;

        newl = m->ciphone[(int) lc].filler ? m->sil : lc;
        newr = m->ciphone[(int) rc].filler ? m->sil : rc;
        if ((newl == lc) && (newr == rc))
            return -1;

        return (mdef_phone_id(m, ci, newl, newr, wpos));
    }

    return (rcptr->pid);
}

int
mdef_is_ciphone(mdef_t * m, int p)
{
    assert(m);
    assert((p >= 0) && (p < m->n_phone));

    return ((p < m->n_ciphone) ? 1 : 0);
}

int
mdef_is_cisenone(mdef_t * m, int s)
{
    assert(m);
    if (s >= m->n_sen) {
        return 0;
    }
    assert(s >= 0);
    return ((s == m->cd2cisen[s]) ? 1 : 0);
}


/* Parse tmat and state->senone mappings for phone p and fill in structure */
static void
parse_tmat_senmap(mdef_t * m, char *line, int32 off, int p)
{
    int32 wlen, n, s;
    __BIGSTACKVARIABLE__ char word[1024], *lp;

    lp = line + off;

    /* Read transition matrix id */
    if ((sscanf(lp, "%d%n", &n, &wlen) != 1) || (n < 0))
        E_FATAL("Missing or bad transition matrix id: %s\n", line);
    m->phone[p].tmat = n;
    if (m->n_tmat <= n)
        E_FATAL("tmat-id(%d) > #tmat in header(%d): %s\n", n, m->n_tmat,
                line);
    lp += wlen;

    /* Read senone mappings for each emitting state */
    for (n = 0; n < m->n_emit_state; n++) {
        if ((sscanf(lp, "%d%n", &s, &wlen) != 1) || (s < 0))
            E_FATAL("Missing or bad state[%d]->senone mapping: %s\n", n,
                    line);

        if ((p < m->n_ciphone) && (m->n_ci_sen <= s))
            E_FATAL("CI-senone-id(%d) > #CI-senones(%d): %s\n", s,
                    m->n_ci_sen, line);
        if (m->n_sen <= s)
            E_FATAL("Senone-id(%d) > #senones(%d): %s\n", s, m->n_sen,
                    line);

        m->sseq[p][n] = s;
        lp += wlen;
    }

    /* Check for the last non-emitting state N */
    if ((sscanf(lp, "%s%n", word, &wlen) != 1) || (strcmp(word, "N") != 0))
        E_FATAL("Missing non-emitting state spec: %s\n", line);
    lp += wlen;

    /* Check for end of line */
    if (sscanf(lp, "%s%n", word, &wlen) == 1)
        E_FATAL("Non-empty beyond non-emitting final state: %s\n", line);
}


static void
parse_base_line(mdef_t * m, char *line, int p)
{
    int32 wlen, n;
    __BIGSTACKVARIABLE__ char word[1024], *lp;
    int ci;

    lp = line;

    /* Read base phone name */
    if (sscanf(lp, "%s%n", word, &wlen) != 1)
        E_FATAL("Missing base phone name: %s\n", line);
    lp += wlen;

    /* Make sure it's not a duplicate */
    ci = mdef_ciphone_id(m, word);
    if (ci >= 0)
        E_FATAL("Duplicate base phone: %s\n", line);

    /* Add ciphone to ciphone table with id p */
    ciphone_add(m, word, p);
    ci = (int) p;

    /* Read and skip "-" for lc, rc, wpos */
    for (n = 0; n < 3; n++) {
        if ((sscanf(lp, "%s%n", word, &wlen) != 1)
            || (strcmp(word, "-") != 0))
            E_FATAL("Bad context info for base phone: %s\n", line);
        lp += wlen;
    }

    /* Read filler attribute, if present */
    if (sscanf(lp, "%s%n", word, &wlen) != 1)
        E_FATAL("Missing filler atribute field: %s\n", line);
    lp += wlen;
    if (strcmp(word, "filler") == 0)
        m->ciphone[(int) ci].filler = 1;
    else if (strcmp(word, "n/a") == 0)
        m->ciphone[(int) ci].filler = 0;
    else
        E_FATAL("Bad filler attribute field: %s\n", line);

    triphone_add(m, ci, -1, -1, WORD_POSN_UNDEFINED, p);

    /* Parse remainder of line: transition matrix and state->senone mappings */
    parse_tmat_senmap(m, line, lp - line, p);
}


static void
parse_tri_line(mdef_t * m, char *line, int p)
{
    int32 wlen;
    __BIGSTACKVARIABLE__ char word[1024], *lp;
    int ci, lc, rc;
    word_posn_t wpos = WORD_POSN_BEGIN;

    lp = line;

    /* Read base phone name */
    if (sscanf(lp, "%s%n", word, &wlen) != 1)
        E_FATAL("Missing base phone name: %s\n", line);
    lp += wlen;

    ci = mdef_ciphone_id(m, word);
    if (ci < 0)
        E_FATAL("Unknown base phone: %s\n", line);

    /* Read lc */
    if (sscanf(lp, "%s%n", word, &wlen) != 1)
        E_FATAL("Missing left context: %s\n", line);
    lp += wlen;
    lc = mdef_ciphone_id(m, word);
    if (lc < 0)
        E_FATAL("Unknown left context: %s\n", line);

    /* Read rc */
    if (sscanf(lp, "%s%n", word, &wlen) != 1)
        E_FATAL("Missing right context: %s\n", line);
    lp += wlen;
    rc = mdef_ciphone_id(m, word);
    if (rc < 0)
        E_FATAL("Unknown right  context: %s\n", line);

    /* Read tripone word-position within word */
    if ((sscanf(lp, "%s%n", word, &wlen) != 1) || (word[1] != '\0'))
        E_FATAL("Missing or bad word-position spec: %s\n", line);
    lp += wlen;
    switch (word[0]) {
    case 'b':
        wpos = WORD_POSN_BEGIN;
        break;
    case 'e':
        wpos = WORD_POSN_END;
        break;
    case 's':
        wpos = WORD_POSN_SINGLE;
        break;
    case 'i':
        wpos = WORD_POSN_INTERNAL;
        break;
    default:
        E_FATAL("Bad word-position spec: %s\n", line);
    }

    /* Read filler attribute, if present.  Must match base phone attribute */
    if (sscanf(lp, "%s%n", word, &wlen) != 1)
        E_FATAL("Missing filler attribute field: %s\n", line);
    lp += wlen;
    if (((strcmp(word, "filler") == 0) && (m->ciphone[(int) ci].filler)) ||
        ((strcmp(word, "n/a") == 0) && (!m->ciphone[(int) ci].filler))) {
        /* Everything is fine */
    }
    else
        E_FATAL("Bad filler attribute field: %s\n", line);

    triphone_add(m, ci, lc, rc, wpos, p);

    /* Parse remainder of line: transition matrix and state->senone mappings */
    parse_tmat_senmap(m, line, lp - line, p);
}


static void
sseq_compress(mdef_t * m)
{
    hash_table_t *h;
    uint16 **sseq;
    int32 n_sseq;
    int32 p, j, k;
    glist_t g;
    gnode_t *gn;
    hash_entry_t *he;

    k = m->n_emit_state * sizeof(int16);

    h = hash_table_new(m->n_phone, HASH_CASE_YES);
    n_sseq = 0;

    /* Identify unique senone-sequence IDs.  BUG: tmat-id not being considered!! */
    for (p = 0; p < m->n_phone; p++) {
        /* Add senone sequence to hash table */
	if (n_sseq
            == (j = hash_table_enter_bkey_int32(h, (char *)m->sseq[p], k, n_sseq)))
            n_sseq++;

        m->phone[p].ssid = j;
    }

    /* Generate compacted sseq table */
    sseq = ckd_calloc_2d(n_sseq, m->n_emit_state, sizeof(**sseq)); /* freed in mdef_free() */

    g = hash_table_tolist(h, &j);
    assert(j == n_sseq);

    for (gn = g; gn; gn = gnode_next(gn)) {
        he = (hash_entry_t *) gnode_ptr(gn);
        j = (long)hash_entry_val(he);
        memcpy(sseq[j], hash_entry_key(he), k);
    }
    glist_free(g);

    /* Free the old, temporary senone sequence table, replace with compacted one */
    ckd_free_2d(m->sseq);
    m->sseq = sseq;
    m->n_sseq = n_sseq;

    hash_table_free(h);
}


static int32
noncomment_line(char *line, int32 size, FILE * fp)
{
    while (fgets(line, size, fp) != NULL) {
        if (line[0] != '#')
            return 0;
    }
    return -1;
}


/*
 * Initialize phones (ci and triphones) and state->senone mappings from .mdef file.
 */
mdef_t *
mdef_init(char *mdeffile, int32 breport)
{
    FILE *fp;
    int32 n_ci, n_tri, n_map, n;
    __BIGSTACKVARIABLE__ char tag[1024], buf[1024];
    uint16 **senmap;
    int p;
    int32 s, ci, cd;
    mdef_t *m;

    if (!mdeffile)
        E_FATAL("No mdef-file\n");

    if (breport)
        E_INFO("Reading model definition: %s\n", mdeffile);

    m = (mdef_t *) ckd_calloc(1, sizeof(mdef_t));       /* freed in mdef_free */

    if ((fp = fopen(mdeffile, "r")) == NULL)
        E_FATAL_SYSTEM("Failed to open mdef file '%s' for reading", mdeffile);

    if (noncomment_line(buf, sizeof(buf), fp) < 0)
        E_FATAL("Empty file: %s\n", mdeffile);

    if (strncmp(buf, "BMDF", 4) == 0 || strncmp(buf, "FDMB", 4) == 0) {
        E_INFO
            ("Found byte-order mark %.4s, assuming this is a binary mdef file\n",
             buf);
        fclose(fp);
        ckd_free(m);
        return NULL;
    }
    if (strncmp(buf, MODEL_DEF_VERSION, strlen(MODEL_DEF_VERSION)) != 0)
        E_FATAL("Version error: Expecing %s, but read %s\n",
                MODEL_DEF_VERSION, buf);

    /* Read #base phones, #triphones, #senone mappings defined in header */
    n_ci = -1;
    n_tri = -1;
    n_map = -1;
    m->n_ci_sen = -1;
    m->n_sen = -1;
    m->n_tmat = -1;
    do {
        if (noncomment_line(buf, sizeof(buf), fp) < 0)
            E_FATAL("Incomplete header\n");

        if ((sscanf(buf, "%d %s", &n, tag) != 2) || (n < 0))
            E_FATAL("Error in header: %s\n", buf);

        if (strcmp(tag, "n_base") == 0)
            n_ci = n;
        else if (strcmp(tag, "n_tri") == 0)
            n_tri = n;
        else if (strcmp(tag, "n_state_map") == 0)
            n_map = n;
        else if (strcmp(tag, "n_tied_ci_state") == 0)
            m->n_ci_sen = n;
        else if (strcmp(tag, "n_tied_state") == 0)
            m->n_sen = n;
        else if (strcmp(tag, "n_tied_tmat") == 0)
            m->n_tmat = n;
        else
            E_FATAL("Unknown header line: %s\n", buf);
    } while ((n_ci < 0) || (n_tri < 0) || (n_map < 0) ||
             (m->n_ci_sen < 0) || (m->n_sen < 0) || (m->n_tmat < 0));

    if ((n_ci == 0) || (m->n_ci_sen == 0) || (m->n_tmat == 0)
        || (m->n_ci_sen > m->n_sen))
        E_FATAL("%s: Error in header\n", mdeffile);

    /* Check typesize limits */
    if (n_ci >= MAX_INT16)
        E_FATAL("%s: #CI phones (%d) exceeds limit (%d)\n", mdeffile, n_ci,
                MAX_INT16);
    if (n_ci + n_tri >= MAX_INT32) /* Comparison is always false... */
        E_FATAL("%s: #Phones (%d) exceeds limit (%d)\n", mdeffile,
                n_ci + n_tri, MAX_INT32);
    if (m->n_sen >= MAX_INT16)
        E_FATAL("%s: #senones (%d) exceeds limit (%d)\n", mdeffile,
                m->n_sen, MAX_INT16);
    if (m->n_tmat >= MAX_INT32) /* Comparison is always false... */
        E_FATAL("%s: #tmats (%d) exceeds limit (%d)\n", mdeffile,
                m->n_tmat, MAX_INT32);

    m->n_emit_state = (n_map / (n_ci + n_tri)) - 1;
    if ((m->n_emit_state + 1) * (n_ci + n_tri) != n_map)
        E_FATAL
            ("Header error: n_state_map not a multiple of n_ci*n_tri\n");

    /* Initialize ciphone info */
    m->n_ciphone = n_ci;
    m->ciphone_ht = hash_table_new(n_ci, HASH_CASE_YES);  /* With case-insensitive string names *//* freed in mdef_free */
    m->ciphone = (ciphone_t *) ckd_calloc(n_ci, sizeof(ciphone_t));     /* freed in mdef_free */

    /* Initialize phones info (ciphones + triphones) */
    m->n_phone = n_ci + n_tri;
    m->phone = (phone_t *) ckd_calloc(m->n_phone, sizeof(phone_t));     /* freed in mdef_free */

    /* Allocate space for state->senone map for each phone */
    senmap = ckd_calloc_2d(m->n_phone, m->n_emit_state, sizeof(**senmap));      /* freed in mdef_free */
    m->sseq = senmap;           /* TEMPORARY; until it is compressed into just the unique ones */

    /* Allocate initial space for <ci,lc,rc,wpos> -> pid mapping */
    m->wpos_ci_lclist = (ph_lc_t ***) ckd_calloc_2d(N_WORD_POSN, m->n_ciphone, sizeof(ph_lc_t *));      /* freed in mdef_free */

    /*
     * Read base phones and triphones.  They'll simply be assigned a running sequence
     * number as their "phone-id".  If the phone-id < n_ci, it's a ciphone.
     */

    /* Read base phones */
    for (p = 0; p < n_ci; p++) {
        if (noncomment_line(buf, sizeof(buf), fp) < 0)
            E_FATAL("Premature EOF reading CIphone %d\n", p);
        parse_base_line(m, buf, p);
    }
    m->sil = mdef_ciphone_id(m, S3_SILENCE_CIPHONE);

    /* Read triphones, if any */
    for (; p < m->n_phone; p++) {
        if (noncomment_line(buf, sizeof(buf), fp) < 0)
            E_FATAL("Premature EOF reading phone %d\n", p);
        parse_tri_line(m, buf, p);
    }

    if (noncomment_line(buf, sizeof(buf), fp) >= 0)
        E_ERROR("Non-empty file beyond expected #phones (%d)\n",
                m->n_phone);

    /* Build CD senones to CI senones map */
    if (m->n_ciphone * m->n_emit_state != m->n_ci_sen)
        E_FATAL
            ("#CI-senones(%d) != #CI-phone(%d) x #emitting-states(%d)\n",
             m->n_ci_sen, m->n_ciphone, m->n_emit_state);
    m->cd2cisen = (int16 *) ckd_calloc(m->n_sen, sizeof(*m->cd2cisen)); /* freed in mdef_free */

    m->sen2cimap = (int16 *) ckd_calloc(m->n_sen, sizeof(*m->sen2cimap)); /* freed in mdef_free */

    for (s = 0; s < m->n_sen; s++)
        m->sen2cimap[s] = -1;
    for (s = 0; s < m->n_ci_sen; s++) { /* CI senones */
        m->cd2cisen[s] = s;
        m->sen2cimap[s] = s / m->n_emit_state;
    }
    for (p = n_ci; p < m->n_phone; p++) {       /* CD senones */
        for (s = 0; s < m->n_emit_state; s++) {
            cd = m->sseq[p][s];
            ci = m->sseq[m->phone[p].ci][s];
            m->cd2cisen[cd] = ci;
            m->sen2cimap[cd] = m->phone[p].ci;
        }
    }

    sseq_compress(m);
    fclose(fp);

    return m;
}

void
mdef_report(mdef_t * m)
{
    E_INFO_NOFN("Initialization of mdef_t, report:\n");
    E_INFO_NOFN
        ("%d CI-phone, %d CD-phone, %d emitstate/phone, %d CI-sen, %d Sen, %d Sen-Seq\n",
         m->n_ciphone, m->n_phone - m->n_ciphone, m->n_emit_state,
         m->n_ci_sen, m->n_sen, m->n_sseq);
    E_INFO_NOFN("\n");

}

/* RAH 4.23.01, Need to step down the ->next list to see if there are
   any more things to free
 */



/* RAH 4.19.01, Attempt to free memory that was allocated within this module
   I have not verified that all the memory has been freed. I've taken only a 
   reasonable effort for now.
   RAH 4.24.01 - verified that all memory is released.
 */
void
mdef_free_recursive_lc(ph_lc_t * lc)
{
    if (lc == NULL)
        return;

    if (lc->rclist)
        mdef_free_recursive_rc(lc->rclist);

    if (lc->next)
        mdef_free_recursive_lc(lc->next);

    ckd_free((void *) lc);
}

void
mdef_free_recursive_rc(ph_rc_t * rc)
{
    if (rc == NULL)
        return;

    if (rc->next)
        mdef_free_recursive_rc(rc->next);

    ckd_free((void *) rc);
}


/* RAH, Free memory that was allocated in mdef_init 
   Rational purify shows that no leaks exist
 */

void
mdef_free(mdef_t * m)
{
    int i, j;

    if (m) {
        if (m->sen2cimap)
            ckd_free((void *) m->sen2cimap);
        if (m->cd2cisen)
            ckd_free((void *) m->cd2cisen);

        /* RAH, go down the ->next list and delete all the pieces */
        for (i = 0; i < N_WORD_POSN; i++)
            for (j = 0; j < m->n_ciphone; j++)
                if (m->wpos_ci_lclist[i][j]) {
                    mdef_free_recursive_lc(m->wpos_ci_lclist[i][j]->next);
                    mdef_free_recursive_rc(m->wpos_ci_lclist[i][j]->
                                           rclist);
                }

        for (i = 0; i < N_WORD_POSN; i++)
            for (j = 0; j < m->n_ciphone; j++)
                if (m->wpos_ci_lclist[i][j])
                    ckd_free((void *) m->wpos_ci_lclist[i][j]);


        if (m->wpos_ci_lclist)
            ckd_free_2d((void *) m->wpos_ci_lclist);
        if (m->sseq)
            ckd_free_2d((void *) m->sseq);
        /* Free phone context */
        if (m->phone)
            ckd_free((void *) m->phone);
        if (m->ciphone_ht)
            hash_table_free(m->ciphone_ht);

        for (i = 0; i < m->n_ciphone; i++) {
            if (m->ciphone[i].name)
                ckd_free((void *) m->ciphone[i].name);
        }


        if (m->ciphone)
            ckd_free((void *) m->ciphone);

        ckd_free((void *) m);
    }
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * mdef.h -- HMM model definition: base (CI) phones and triphones
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1999 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 */


#ifndef __MDEF_H__
#define __MDEF_H__


/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/hash_table.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

/** \file mdef.h
 * \brief Model definition 
 */

/** \enum word_posn_t
 * \brief Union of different type of word position
 */

typedef enum {
    WORD_POSN_INTERNAL = 0,	/**< Internal phone of word */
    WORD_POSN_BEGIN = 1,	/**< Beginning phone of word */
    WORD_POSN_END = 2,		/**< Ending phone of word */
    WORD_POSN_SINGLE = 3,	/**< Single phone word (i.e. begin & end) */
    WORD_POSN_UNDEFINED = 4	/**< Undefined value, used for initial conditions, etc */
} word_posn_t;
#define N_WORD_POSN	4	/**< total # of word positions (excluding undefined) */
#define WPOS_NAME	"ibesu"	/**< Printable code for each word position above */
#define S3_SILENCE_CIPHONE "SIL" /**< Hard-coded silence CI phone name */

/**
   \struct ciphone_t
   \brief CI phone information 
*/
typedef struct {
    char *name;                 /**< The name of the CI phone */
    int32 filler;		/**< Whether a filler phone; if so, can be substituted by
				   silence phone in left or right context position */
} ciphone_t;

/**
 * \struct phone_t
 * \brief Triphone information, including base phones as a subset.  For the latter, lc, rc and wpos are non-existent.
 */
typedef struct {
    int32 ssid;			/**< State sequence (or senone sequence) ID, considering the
				   n_emit_state senone-ids are a unit.  The senone sequences
				   themselves are in a separate table */
    int32 tmat;			/**< Transition matrix id */
    int16 ci, lc, rc;		/**< Base, left, right context ciphones */
    word_posn_t wpos;		/**< Word position */
    
} phone_t;

/**
 * \struct ph_rc_t
 * \brief Structures needed for mapping <ci,lc,rc,wpos> into pid.  (See mdef_t.wpos_ci_lclist below.)  (lc = left context; rc = right context.)
 * NOTE: Both ph_rc_t and ph_lc_t FOR INTERNAL USE ONLY.
 */
typedef struct ph_rc_s {
    int16 rc;			/**< Specific rc for a parent <wpos,ci,lc> */
    int32 pid;			/**< Triphone id for above rc instance */
    struct ph_rc_s *next;	/**< Next rc entry for same parent <wpos,ci,lc> */
} ph_rc_t;

/**
 * \struct ph_lc_t
 * \brief Structures for storing the left context. 
 */

typedef struct ph_lc_s {
    int16 lc;			/**< Specific lc for a parent <wpos,ci> */
    ph_rc_t *rclist;		/**< rc list for above lc instance */
    struct ph_lc_s *next;	/**< Next lc entry for same parent <wpos,ci> */
} ph_lc_t;


/** The main model definition structure */
/**
   \struct mdef_t
   \brief strcture for storing the model definition. 
*/
typedef struct {
    int32 n_ciphone;		/**< number basephones actually present */
    int32 n_phone;		/**< number basephones + number triphones actually present */
    int32 n_emit_state;		/**< number emitting states per phone */
    int32 n_ci_sen;		/**< number CI senones; these are the first */
    int32 n_sen;		/**< number senones (CI+CD) */
    int32 n_tmat;		/**< number transition matrices */
    
    hash_table_t *ciphone_ht;	/**< Hash table for mapping ciphone strings to ids */
    ciphone_t *ciphone;		/**< CI-phone information for all ciphones */
    phone_t *phone;		/**< Information for all ciphones and triphones */
    uint16 **sseq;		/**< Unique state (or senone) sequences in this model, shared
                                   among all phones/triphones */
    int32 n_sseq;		/**< No. of unique senone sequences in this model */
    
    int16 *cd2cisen;		/**< Parent CI-senone id for each senone; the first
				   n_ci_sen are identity mappings; the CD-senones are
				   contiguous for each parent CI-phone */
    int16 *sen2cimap;		/**< Parent CI-phone for each senone (CI or CD) */
    
    int16 sil;			/**< SILENCE_CIPHONE id */
    
    ph_lc_t ***wpos_ci_lclist;	/**< wpos_ci_lclist[wpos][ci] = list of lc for <wpos,ci>.
                                   wpos_ci_lclist[wpos][ci][lc].rclist = list of rc for
                                   <wpos,ci,lc>.  Only entries for the known triphones
                                   are created to conserve space.
                                   (NOTE: FOR INTERNAL USE ONLY.) */
} mdef_t;

/** Access macros; not meant for arbitrary use */
#define mdef_is_fillerphone(m,p)	((m)->ciphone[p].filler)
#define mdef_n_ciphone(m)		((m)->n_ciphone)
#define mdef_n_phone(m)			((m)->n_phone)
#define mdef_n_sseq(m)			((m)->n_sseq)
#define mdef_n_emit_state(m)		((m)->n_emit_state)
#define mdef_n_sen(m)			((m)->n_sen)
#define mdef_n_tmat(m)			((m)->n_tmat)
#define mdef_pid2ssid(m,p)		((m)->phone[p].ssid)
#define mdef_pid2tmatid(m,p)		((m)->phone[p].tmat)
#define mdef_silphone(m)		((m)->sil)
#define mdef_sen2cimap(m)		((m)->sen2cimap)
#define mdef_sseq2sen(m,ss,pos)		((m)->sseq[ss][pos])
#define mdef_pid2ci(m,p)		((m)->phone[p].ci)
#define mdef_cd2cisen(m)		((m)->cd2cisen)

/**
 * Initialize the phone structure from the given model definition file.
 * It should be treated as a READ-ONLY structure.
 * @return pointer to the phone structure created.
 */
mdef_t *mdef_init (char *mdeffile, /**< In: Model definition file */
		   int breport     /**< In: whether to report the progress or not */
    );


/** 
    Get the ciphone id given a string name
    @return ciphone id for the given ciphone string name 
*/
int mdef_ciphone_id(mdef_t *m,		/**< In: Model structure being queried */
                    char *ciphone	/**< In: ciphone for which id wanted */
    );

/** 
    Get the phone string given the ci phone id.
    @return: READ-ONLY ciphone string name for the given ciphone id 
*/
const char *mdef_ciphone_str(mdef_t *m,	/**< In: Model structure being queried */
                             int ci	/**< In: ciphone id for which name wanted */
    );

/** 
    Decide whether the phone is ci phone.
    @return 1 if given triphone argument is a ciphone, 0 if not, -1 if error 
*/
int mdef_is_ciphone (mdef_t *m,		/**< In: Model structure being queried */
                     int p		/**< In: triphone id being queried */
    );

/**
   Decide whether the senone is a senone for a ci phone, or a ci senone
   @return 1 if a given senone is a ci senone
*/  
int mdef_is_cisenone(mdef_t *m,               /**< In: Model structure being queried */
                     int s		        /**< In: senone id being queried */
    );

/** 
    Decide the phone id given the left, right and base phones. 
    @return: phone id for the given constituents if found, else BAD_S3PID 
*/
int mdef_phone_id (mdef_t *m,		/**< In: Model structure being queried */
                   int b,		/**< In: base ciphone id */
                   int l,		/**< In: left context ciphone id */
                   int r,		/**< In: right context ciphone id */
                   word_posn_t pos	/**< In: Word position */
    );

/**
 * Create a phone string for the given phone (base or triphone) id in the given buf.
 * @return 0 if successful, -1 if error.
 */
int mdef_phone_str(mdef_t *m,		/**< In: Model structure being queried */
                   int pid,		/**< In: phone id being queried */
                   char *buf		/**< Out: On return, buf has the string */
    );

/**
 * Compare the underlying HMMs for two given phones (i.e., compare the two transition
 * matrix IDs and the individual state(senone) IDs).
 * @return 0 iff the HMMs are identical, -1 otherwise.
 */
int mdef_hmm_cmp (mdef_t *m,	/**< In: Model being queried */
                  int p1, 	/**< In: One of the two triphones being compared */
                  int p2	/**< In: One of the two triphones being compared */
    );

/** Report the model definition's parameters */
void mdef_report(mdef_t *m /**<  In: model definition structure */
    );

/** RAH, For freeing memory */
void mdef_free_recursive_lc (ph_lc_t *lc /**< In: A list of left context */
    );
void mdef_free_recursive_rc (ph_rc_t *rc /**< In: A list of right context */
    );

/** Free an mdef_t */
void mdef_free (mdef_t *mdef /**< In : The model definition*/
    );


#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif

#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * gauden.c -- gaussian density module.
 *
 ***********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1996 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 ***********************************************
 *
 * HISTORY
 * $Log$
 * Revision 1.7  2006/02/22  17:09:55  arthchan2003
 * Merged from SPHINX3_5_2_RCI_IRII_BRANCH: 1, Followed Dave's change, keep active to be uint8 instead int8 in gauden_dist_norm.\n 2, Introdued gauden_dump and gauden_dump_ind.  This allows debugging of ms_gauden routine. \n 3, Introduced gauden_free, this fixed some minor memory leaks. \n 4, gauden_init accept an argument precompute to specify whether the distance is pre-computed or not.\n 5, Added license. \n 6, Fixed dox-doc.
 * 
 *
 * Revision 1.5.4.7  2006/01/16 19:45:59  arthchan2003
 * Change the gaussian density dumping routine to a function.
 *
 * Revision 1.5.4.6  2005/10/09 19:51:05  arthchan2003
 * Followed Dave's changed in the trunk.
 *
 * Revision 1.5.4.5  2005/09/25 18:54:20  arthchan2003
 * Added a flag to turn on and off precomputation.
 *
 * Revision 1.6  2005/10/05 00:31:14  dhdfu
 * Make int8 be explicitly signed (signedness of 'char' is
 * architecture-dependent).  Then make a bunch of things use uint8 where
 * signedness is unimportant, because on the architecture where 'char' is
 * unsigned, it is that way for a reason (signed chars are slower).
 *
 * Revision 1.5.4.4  2005/09/07 23:29:07  arthchan2003
 * Added FIXME warning.
 *
 * Revision 1.5.4.3  2005/09/07 23:25:10  arthchan2003
 * 1, Behavior changes of cont_mgau, instead of remove Gaussian with zero variance vector before flooring, now remove Gaussian with zero mean and variance before flooring. Notice that this is not yet synchronize with ms_mgau. 2, Added warning message in multi-stream gaussian distribution.
 *
 * Revision 1.5.4.2  2005/08/03 18:53:44  dhdfu
 * Add memory deallocation functions.  Also move all the initialization
 * of ms_mgau_model_t into ms_mgau_init (duh!), which entails removing it
 * from decode_anytopo and friends.
 *
 * Revision 1.5.4.1  2005/07/20 19:39:01  arthchan2003
 * Added licences in ms_* series of code.
 *
 * Revision 1.5  2005/06/21 18:55:09  arthchan2003
 * 1, Add comments to describe this modules, 2, Fixed doxygen documentation. 3, Added $ keyword.
 *
 * Revision 1.3  2005/03/30 01:22:47  archan
 * Fixed mistakes in last updates. Add
 *
 * 
 * 20-Dec-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Changed gauden_param_read to use the new libio/bio_fread functions.
 * 
 * 26-Sep-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added gauden_mean_reload() for application of MLLR; and correspondingly
 * 		made gauden_param_read allocate memory for parameter only if not
 * 		already allocated.
 * 
 * 09-Sep-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Interleaved two density computations for speed improvement.
 * 
 * 19-Aug-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added compute_dist_all special case for improving speed.
 * 
 * 26-Jan-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added check for underflow and floor insertion in gauden_dist.
 * 
 * 20-Jan-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added active argument to gauden_dist_norm and gauden_dist_norm_global,
 * 		and made the latter a static function.
 * 
 * 07-Nov-95	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Initial version created.
 * 		Very liberally borrowed/adapted from Eric's S3 trainer implementation.
 */

/* System headers. */
#include <assert.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* SphinxBase headers. */
#include <sphinxbase/bio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>

/* Local headesr. */
#include "ms_gauden.h"

#define GAUDEN_PARAM_VERSION	"1.0"

#ifndef M_PI
#define M_PI	3.1415926535897932385e0
#endif

#define WORST_DIST	(int32)(0x80000000)

void
gauden_dump(const gauden_t * g)
{
    int32 c;

    for (c = 0; c < g->n_mgau; c++)
        gauden_dump_ind(g, c);
}


void
gauden_dump_ind(const gauden_t * g, int senidx)
{
    int32 f, d, i;

    for (f = 0; f < g->n_feat; f++) {
        E_INFO("Codebook %d, Feature %d (%dx%d):\n",
               senidx, f, g->n_density, g->featlen[f]);

        for (d = 0; d < g->n_density; d++) {
            printf("m[%3d]", d);
            for (i = 0; i < g->featlen[f]; i++)
		printf(" %7.4f", MFCC2FLOAT(g->mean[senidx][f][d][i]));
            printf("\n");
        }
        printf("\n");

        for (d = 0; d < g->n_density; d++) {
            printf("v[%3d]", d);
            for (i = 0; i < g->featlen[f]; i++)
                printf(" %d", (int)g->var[senidx][f][d][i]);
            printf("\n");
        }
        printf("\n");

        for (d = 0; d < g->n_density; d++)
            printf("d[%3d] %d\n", d, (int)g->det[senidx][f][d]);
    }
    fflush(stderr);
}

static int32
gauden_param_read(float32 ***** out_param,      /* Alloc space iff *out_param == NULL */
                  int32 * out_n_mgau,
                  int32 * out_n_feat,
                  int32 * out_n_density,
                  int32 ** out_veclen, const char *file_name)
{
    char tmp;
    FILE *fp;
    int32 i, j, k, l, n, blk;
    int32 n_mgau;
    int32 n_feat;
    int32 n_density;
    int32 *veclen;
    int32 byteswap, chksum_present;
    float32 ****out;
    float32 *buf;
    char **argname, **argval;
    uint32 chksum;

    E_INFO("Reading mixture gaussian parameter: %s\n", file_name);

    if ((fp = fopen(file_name, "rb")) == NULL)
        E_FATAL_SYSTEM("Failed to open file '%s' for reading", file_name);

    /* Read header, including argument-value info and 32-bit byteorder magic */
    if (bio_readhdr(fp, &argname, &argval, &byteswap) < 0)
        E_FATAL("Failed to read header from file '%s'\n", file_name);

    /* Parse argument-value list */
    chksum_present = 0;
    for (i = 0; argname[i]; i++) {
        if (strcmp(argname[i], "version") == 0) {
            if (strcmp(argval[i], GAUDEN_PARAM_VERSION) != 0)
                E_WARN("Version mismatch(%s): %s, expecting %s\n",
                       file_name, argval[i], GAUDEN_PARAM_VERSION);
        }
        else if (strcmp(argname[i], "chksum0") == 0) {
            chksum_present = 1; /* Ignore the associated value */
        }
    }
    bio_hdrarg_free(argname, argval);
    argname = argval = NULL;

    chksum = 0;

    /* #Codebooks */
    if (bio_fread(&n_mgau, sizeof(int32), 1, fp, byteswap, &chksum) != 1)
        E_FATAL("fread(%s) (#codebooks) failed\n", file_name);
    *out_n_mgau = n_mgau;

    /* #Features/codebook */
    if (bio_fread(&n_feat, sizeof(int32), 1, fp, byteswap, &chksum) != 1)
        E_FATAL("fread(%s) (#features) failed\n", file_name);
    *out_n_feat = n_feat;

    /* #Gaussian densities/feature in each codebook */
    if (bio_fread(&n_density, sizeof(int32), 1, fp, byteswap, &chksum) != 1)
        E_FATAL("fread(%s) (#density/codebook) failed\n", file_name);
    *out_n_density = n_density;

    /* #Dimensions in each feature stream */
    veclen = ckd_calloc(n_feat, sizeof(uint32));
    *out_veclen = veclen;
    if (bio_fread(veclen, sizeof(int32), n_feat, fp, byteswap, &chksum) !=
        n_feat)
        E_FATAL("fread(%s) (feature-lengths) failed\n", file_name);

    /* blk = total vector length of all feature streams */
    for (i = 0, blk = 0; i < n_feat; i++)
        blk += veclen[i];

    /* #Floats to follow; for the ENTIRE SET of CODEBOOKS */
    if (bio_fread(&n, sizeof(int32), 1, fp, byteswap, &chksum) != 1)
        E_FATAL("fread(%s) (total #floats) failed\n", file_name);
    if (n != n_mgau * n_density * blk) {
        E_FATAL
            ("%s: #mfcc_ts(%d) doesn't match dimensions: %d x %d x %d\n",
             file_name, n, n_mgau, n_density, blk);
    }

    /* Allocate memory for mixture gaussian densities if not already allocated */
    if (!(*out_param)) {
        out = (float32 ****) ckd_calloc_3d(n_mgau, n_feat, n_density,
                                         sizeof(float32 *));
        buf = (float32 *) ckd_calloc(n, sizeof(float32));
        for (i = 0, l = 0; i < n_mgau; i++) {
            for (j = 0; j < n_feat; j++) {
                for (k = 0; k < n_density; k++) {
                    out[i][j][k] = &buf[l];
                    l += veclen[j];
                }
            }
        }
    }
    else {
        out = (float32 ****) *out_param;
        buf = out[0][0][0];
    }

    /* Read mixture gaussian densities data */
    if (bio_fread(buf, sizeof(float32), n, fp, byteswap, &chksum) != n)
        E_FATAL("fread(%s) (densitydata) failed\n", file_name);

    if (chksum_present)
        bio_verify_chksum(fp, byteswap, chksum);

    if (fread(&tmp, 1, 1, fp) == 1)
        E_FATAL("More data than expected in %s\n", file_name);

    fclose(fp);

    *out_param = out;

    E_INFO("%d codebook, %d feature, size: \n", n_mgau, n_feat);
    for (i = 0; i < n_feat; i++)
        E_INFO(" %dx%d\n", n_density, veclen[i]);

    return 0;
}

static void
gauden_param_free(mfcc_t **** p)
{
    ckd_free(p[0][0][0]);
    ckd_free_3d(p);
}

/*
 * Some of the gaussian density computation can be carried out in advance:
 * 	log(determinant) calculation,
 * 	1/(2*var) in the exponent,
 * NOTE; The density computation is performed in log domain.
 */
static int32
gauden_dist_precompute(gauden_t * g, logmath_t *lmath, float32 varfloor)
{
    int32 i, m, f, d, flen;
    mfcc_t *meanp;
    mfcc_t *varp;
    mfcc_t *detp;
    int32 floored;

    floored = 0;
    /* Allocate space for determinants */
    g->det = ckd_calloc_3d(g->n_mgau, g->n_feat, g->n_density, sizeof(***g->det));

    for (m = 0; m < g->n_mgau; m++) {
        for (f = 0; f < g->n_feat; f++) {
            flen = g->featlen[f];

            /* Determinants for all variance vectors in g->[m][f] */
            for (d = 0, detp = g->det[m][f]; d < g->n_density; d++, detp++) {
                *detp = 0;
                for (i = 0, varp = g->var[m][f][d], meanp = g->mean[m][f][d];
                     i < flen; i++, varp++, meanp++) {
                    float32 *fvarp = (float32 *)varp;

#ifdef FIXED_POINT
                    float32 *fmp = (float32 *)meanp;
                    *meanp = FLOAT2MFCC(*fmp);
#endif
                    if (*fvarp < varfloor) {
                        *fvarp = varfloor;
                        ++floored;
                    }
                    *detp += (mfcc_t)logmath_log(lmath,
                                                 1.0 / sqrt(*fvarp * 2.0 * M_PI));
                    /* Precompute this part of the exponential */
                    *varp = (mfcc_t)logmath_ln_to_log(lmath,
                                                      (1.0 / (*fvarp * 2.0)));
                }
            }
        }
    }

    E_INFO("%d variance values floored\n", floored);

    return 0;
}


gauden_t *
gauden_init(char const *meanfile, char const *varfile, float32 varfloor, logmath_t *lmath)
{
    int32 i, m, f, d, *flen;
    float32 ****fgau;
    gauden_t *g;

    assert(meanfile != NULL);
    assert(varfile != NULL);
    assert(varfloor > 0.0);

    g = (gauden_t *) ckd_calloc(1, sizeof(gauden_t));
    g->lmath = lmath;

    /* Read means and (diagonal) variances for all mixture gaussians */
    fgau = NULL;
    gauden_param_read(&fgau, &g->n_mgau, &g->n_feat, &g->n_density,
                      &g->featlen, meanfile);
    g->mean = (mfcc_t ****)fgau;
    fgau = NULL;
    gauden_param_read(&fgau, &m, &f, &d, &flen, varfile);
    g->var = (mfcc_t ****)fgau;

    /* Verify mean and variance parameter dimensions */
    if ((m != g->n_mgau) || (f != g->n_feat) || (d != g->n_density))
        E_FATAL
            ("Mixture-gaussians dimensions for means and variances differ\n");
    for (i = 0; i < g->n_feat; i++)
        if (g->featlen[i] != flen[i])
            E_FATAL("Feature lengths for means and variances differ\n");
    ckd_free(flen);

    /* Floor variances and precompute variance determinants */
    gauden_dist_precompute(g, lmath, varfloor);

    return g;
}

void
gauden_free(gauden_t * g)
{
    if (g == NULL)
        return;
    if (g->mean)
        gauden_param_free(g->mean);
    if (g->var)
        gauden_param_free(g->var);
    if (g->det)
        ckd_free_3d(g->det);
    if (g->featlen)
        ckd_free(g->featlen);
    ckd_free(g);
}

/* See compute_dist below */
static int32
compute_dist_all(gauden_dist_t * out_dist, mfcc_t* obs, int32 featlen,
                 mfcc_t ** mean, mfcc_t ** var, mfcc_t * det,
                 int32 n_density)
{
    int32 i, d;

    for (d = 0; d < n_density; ++d) {
        mfcc_t *m;
        mfcc_t *v;
        mfcc_t dval;

        m = mean[d];
        v = var[d];
        dval = det[d];

        for (i = 0; i < featlen; i++) {
            mfcc_t diff;
#ifdef FIXED_POINT
            /* Have to check for underflows here. */
            mfcc_t pdval = dval;
            diff = obs[i] - m[i];
            dval -= MFCCMUL(MFCCMUL(diff, diff), v[i]);
            if (dval > pdval) {
                dval = WORST_SCORE;
                break;
            }
#else
            diff = obs[i] - m[i];
            /* The compiler really likes this to be a single
             * expression, for whatever reason. */
            dval -= diff * diff * v[i];
#endif
        }

        out_dist[d].dist = dval;
        out_dist[d].id = d;
    }

    return 0;
}


/*
 * Compute the top-N closest gaussians from the chosen set (mgau,feat)
 * for the given input observation vector.
 */
static int32
compute_dist(gauden_dist_t * out_dist, int32 n_top,
             mfcc_t * obs, int32 featlen,
             mfcc_t ** mean, mfcc_t ** var, mfcc_t * det,
             int32 n_density)
{
    int32 i, j, d;
    gauden_dist_t *worst;

    /* Special case optimization when n_density <= n_top */
    if (n_top >= n_density)
        return (compute_dist_all
                (out_dist, obs, featlen, mean, var, det, n_density));

    for (i = 0; i < n_top; i++)
        out_dist[i].dist = WORST_DIST;
    worst = &(out_dist[n_top - 1]);

    for (d = 0; d < n_density; d++) {
        mfcc_t *m;
        mfcc_t *v;
        mfcc_t dval;

        m = mean[d];
        v = var[d];
        dval = det[d];

        for (i = 0; (i < featlen) && (dval >= worst->dist); i++) {
            mfcc_t diff;
#ifdef FIXED_POINT
            /* Have to check for underflows here. */
            mfcc_t pdval = dval;
            diff = obs[i] - m[i];
            dval -= MFCCMUL(MFCCMUL(diff, diff), v[i]);
            if (dval > pdval) {
                dval = WORST_SCORE;
                break;
            }
#else
            diff = obs[i] - m[i];
            /* The compiler really likes this to be a single
             * expression, for whatever reason. */
            dval -= diff * diff * v[i];
#endif
        }

        if ((i < featlen) || (dval < worst->dist))     /* Codeword d worse than worst */
            continue;

        /* Codeword d at least as good as worst so far; insert in the ordered list */
        for (i = 0; (i < n_top) && (dval < out_dist[i].dist); i++);
        assert(i < n_top);
        for (j = n_top - 1; j > i; --j)
            out_dist[j] = out_dist[j - 1];
        out_dist[i].dist = dval;
        out_dist[i].id = d;
    }

    return 0;
}


/*
 * Compute distances of the input observation from the top N codewords in the given
 * codebook (g->{mean,var}[mgau]).  The input observation, obs, includes vectors for
 * all features in the codebook.
 */
int32
gauden_dist(gauden_t * g,
            int mgau, int32 n_top, mfcc_t** obs, gauden_dist_t ** out_dist)
{
    int32 f;

    assert((n_top > 0) && (n_top <= g->n_density));

    for (f = 0; f < g->n_feat; f++) {
        compute_dist(out_dist[f], n_top,
                     obs[f], g->featlen[f],
                     g->mean[mgau][f], g->var[mgau][f], g->det[mgau][f],
                     g->n_density);
        E_DEBUG(3, ("Top CW(%d,%d) = %d %d\n", mgau, f, out_dist[f][0].id,
                    (int)out_dist[f][0].dist >> SENSCR_SHIFT));
    }

    return 0;
}

int32
gauden_mllr_transform(gauden_t *g, ps_mllr_t *mllr, cmd_ln_t *config)
{
    int32 i, m, f, d, *flen;
    float32 ****fgau;

    /* Free data if already here */
    if (g->mean)
        gauden_param_free(g->mean);
    if (g->var)
        gauden_param_free(g->var);
    if (g->det)
        ckd_free_3d(g->det);
    if (g->featlen)
        ckd_free(g->featlen);
    g->mean = NULL;
    g->var = NULL;
    g->det = NULL;
    g->featlen = NULL;

    /* Reload means and variances (un-precomputed). */
    fgau = NULL;
    gauden_param_read(&fgau, &g->n_mgau, &g->n_feat, &g->n_density,
                      &g->featlen, cmd_ln_str_r(config, "-mean"));
    g->mean = (mfcc_t ****)fgau;
    fgau = NULL;
    gauden_param_read(&fgau, &m, &f, &d, &flen, cmd_ln_str_r(config, "-var"));
    g->var = (mfcc_t ****)fgau;

    /* Verify mean and variance parameter dimensions */
    if ((m != g->n_mgau) || (f != g->n_feat) || (d != g->n_density))
        E_FATAL
            ("Mixture-gaussians dimensions for means and variances differ\n");
    for (i = 0; i < g->n_feat; i++)
        if (g->featlen[i] != flen[i])
            E_FATAL("Feature lengths for means and variances differ\n");
    ckd_free(flen);

    /* Transform codebook for each stream s */
    for (i = 0; i < g->n_mgau; ++i) {
        for (f = 0; f < g->n_feat; ++f) {
            float64 *temp;
            temp = (float64 *) ckd_calloc(g->featlen[f], sizeof(float64));
            /* Transform each density d in selected codebook */
            for (d = 0; d < g->n_density; d++) {
                int l;
                for (l = 0; l < g->featlen[f]; l++) {
                    temp[l] = 0.0;
                    for (m = 0; m < g->featlen[f]; m++) {
                        /* FIXME: For now, only one class, hence the zeros below. */
                        temp[l] += mllr->A[f][0][l][m] * g->mean[i][f][d][m];
                    }
                    temp[l] += mllr->b[f][0][l];
                }

                for (l = 0; l < g->featlen[f]; l++) {
                    g->mean[i][f][d][l] = (float32) temp[l];
                    g->var[i][f][d][l] *= mllr->h[f][0][l];
                }
            }
            ckd_free(temp);
        }
    }

    /* Re-precompute (if we aren't adapting variances this isn't
     * actually necessary...) */
    gauden_dist_precompute(g, g->lmath, cmd_ln_float32_r(config, "-varfloor"));
    return 0;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

#ifndef _LIBFBS_GAUDEN_H_
#define _LIBFBS_GAUDEN_H_

/** \file ms_gauden.h
 * \brief (Sphinx 3.0 specific) Gaussian density module.
 *
 * Gaussian density distribution implementation. There are two major
 * difference bettwen ms_gauden and cont_mgau. One is the fact that
 * ms_gauden only take cares of the Gaussian computation part where
 * cont_mgau actually take care of senone computation as well. The
 * other is the fact that ms_gauden is a multi-stream implementation
 * of GMM computation.
 *
 */

/* SphinxBase headers. */
#include <sphinxbase/feat.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/cmd_ln.h>

/* Local headers. */
#include "vector.h"
#include "pocketsphinx_internal.h"
#include "hmm.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

/**
 * \struct gauden_dist_t
 * \brief Structure to store distance (density) values for a given input observation wrt density values in some given codebook.
 */
typedef struct {
    int32 id;		/**< Index of codeword (gaussian density) */
    mfcc_t dist;		/**< Density value for input observation wrt above codeword;
                           NOTE: result in logs3 domain, but var_t used for speed */

} gauden_dist_t;

/**
 * \struct gauden_t
 * \brief Multivariate gaussian mixture density parameters
 */
typedef struct {
    mfcc_t ****mean;	/**< mean[codebook][feature][codeword] vector */
    mfcc_t ****var;	/**< like mean; diagonal covariance vector only */
    mfcc_t ***det;	/**< log(determinant) for each variance vector;
			   actually, log(sqrt(2*pi*det)) */
    logmath_t *lmath;   /**< log math computation */
    int32 n_mgau;	/**< Number codebooks */
    int32 n_feat;	/**< Number feature streams in each codebook */
    int32 n_density;	/**< Number gaussian densities in each codebook-feature stream */
    int32 *featlen;	/**< feature length for each feature */
} gauden_t;


/**
 * Read mixture gaussian codebooks from the given files.  Allocate memory space needed
 * for them.  Apply the specified variance floor value.
 * Return value: ptr to the model created; NULL if error.
 * (See Sphinx3 model file-format documentation.)
 */
gauden_t *
gauden_init (char const *meanfile,/**< Input: File containing means of mixture gaussians */
	     char const *varfile,/**< Input: File containing variances of mixture gaussians */
	     float32 varfloor,	/**< Input: Floor value to be applied to variances */
             logmath_t *lmath
    );

/** Release memory allocated by gauden_init. */
void gauden_free(gauden_t *g); /**< In: The gauden_t to free */

/** Transform Gaussians according to an MLLR matrix (or, eventually, more). */
int32 gauden_mllr_transform(gauden_t *s, ps_mllr_t *mllr, cmd_ln_t *config);

/**
 * Compute gaussian density values for the given input observation vector wrt the
 * specified mixture gaussian codebook (which may consist of several feature streams).
 * Density values are left UNnormalized.
 * @return 0 if successful, -1 otherwise.
 */
int32
gauden_dist (gauden_t *g,	/**< In: handle to entire ensemble of codebooks */
	     int mgau,		/**< In: codebook for which density values to be evaluated
				   (g->{mean,var}[mgau]) */
	     int n_top,		/**< In: Number top densities to be evaluated */
	     mfcc_t **obs,	/**< In: Observation vector; obs[f] = for feature f */
	     gauden_dist_t **out_dist
	     /**< Out: n_top best codewords and density values,
		in worsening order, for each feature stream.
		out_dist[f][i] = i-th best density for feature f.
		Caller must allocate memory for this output */
    );

/**
   Dump the definitionn of Gaussian distribution. 
*/
void gauden_dump (const gauden_t *g  /**< In: Gaussian distribution g*/
    );

/**
   Dump the definition of Gaussian distribution of a particular index to the standard output stream
*/
void gauden_dump_ind (const gauden_t *g,  /**< In: Gaussian distribution g*/
		      int senidx          /**< In: The senone index of the Gaussian */
    );

#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif

#endif /* GAUDEN_H */ 
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * ms_mgau.c -- Essentially a wrapper that wrap up gauden and
 * senone. It supports multi-stream. 
 *
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1997 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * HISTORY
 * $Log$
 * Revision 1.2  2006/02/22  16:56:01  arthchan2003
 * Merged from SPHINX3_5_2_RCI_IRII_BRANCH: Added ms_mgau.[ch] into the trunk. It is a wrapper of ms_gauden and ms_senone
 * 
 * Revision 1.1.2.4  2005/09/25 18:55:19  arthchan2003
 * Added a flag to turn on and off precomputation.
 *
 * Revision 1.1.2.3  2005/08/03 18:53:44  dhdfu
 * Add memory deallocation functions.  Also move all the initialization
 * of ms_mgau_model_t into ms_mgau_init (duh!), which entails removing it
 * from decode_anytopo and friends.
 *
 * Revision 1.1.2.2  2005/08/02 21:05:38  arthchan2003
 * 1, Added dist and mgau_active as intermediate variable for computation. 2, Added ms_cont_mgau_frame_eval, which is a multi stream version of GMM computation mainly s3.0 family of tools. 3, Fixed dox-doc.
 *
 * Revision 1.1.2.1  2005/07/20 19:37:09  arthchan2003
 * Added a multi-stream cont_mgau (ms_mgau) which is a wrapper of both gauden and senone.  Add ms_mgau_init and model_set_mllr.  This allow eliminating 600 lines of code in decode_anytopo/align/allphone.
 *
 *
 *
 */

/* Local headers. */
#include "ms_mgau.h"

static ps_mgaufuncs_t ms_mgau_funcs = {
    "ms",
    ms_cont_mgau_frame_eval, /* frame_eval */
    ms_mgau_mllr_transform,  /* transform */
    ms_mgau_free             /* free */
};

ps_mgau_t *
ms_mgau_init(acmod_t *acmod, logmath_t *lmath, bin_mdef_t *mdef)
{
    /* Codebooks */
    ms_mgau_model_t *msg;
    ps_mgau_t *mg;
    gauden_t *g;
    senone_t *s;
    cmd_ln_t *config;
    int i;

    config = acmod->config;

    msg = (ms_mgau_model_t *) ckd_calloc(1, sizeof(ms_mgau_model_t));
    msg->config = config;
    msg->g = NULL;
    msg->s = NULL;
    
    g = msg->g = gauden_init(cmd_ln_str_r(config, "-mean"),
                             cmd_ln_str_r(config, "-var"),
                             cmd_ln_float32_r(config, "-varfloor"),
                             lmath);

    /* Verify n_feat and veclen, against acmod. */
    if (g->n_feat != feat_dimension1(acmod->fcb)) {
        E_ERROR("Number of streams does not match: %d != %d\n",
                g->n_feat, feat_dimension1(acmod->fcb));
        goto error_out;
    }
    for (i = 0; i < g->n_feat; ++i) {
        if (g->featlen[i] != feat_dimension2(acmod->fcb, i)) {
            E_ERROR("Dimension of stream %d does not match: %d != %d\n", i,
                    g->featlen[i], feat_dimension2(acmod->fcb, i));
            goto error_out;
        }
    }

    s = msg->s = senone_init(msg->g,
                             cmd_ln_str_r(config, "-mixw"),
                             cmd_ln_str_r(config, "-senmgau"),
                             cmd_ln_float32_r(config, "-mixwfloor"),
                             lmath, mdef);

    s->aw = cmd_ln_int32_r(config, "-aw");

    /* Verify senone parameters against gauden parameters */
    if (s->n_feat != g->n_feat)
        E_FATAL("#Feature mismatch: gauden= %d, senone= %d\n", g->n_feat,
                s->n_feat);
    if (s->n_cw != g->n_density)
        E_FATAL("#Densities mismatch: gauden= %d, senone= %d\n",
                g->n_density, s->n_cw);
    if (s->n_gauden > g->n_mgau)
        E_FATAL("Senones need more codebooks (%d) than present (%d)\n",
                s->n_gauden, g->n_mgau);
    if (s->n_gauden < g->n_mgau)
        E_ERROR("Senones use fewer codebooks (%d) than present (%d)\n",
                s->n_gauden, g->n_mgau);

    msg->topn = cmd_ln_int32_r(config, "-topn");
    E_INFO("The value of topn: %d\n", msg->topn);
    if (msg->topn == 0 || msg->topn > msg->g->n_density) {
        E_WARN
            ("-topn argument (%d) invalid or > #density codewords (%d); set to latter\n",
             msg->topn, msg->g->n_density);
        msg->topn = msg->g->n_density;
    }

    msg->dist = (gauden_dist_t ***)
        ckd_calloc_3d(g->n_mgau, g->n_feat, msg->topn,
                      sizeof(gauden_dist_t));
    msg->mgau_active = ckd_calloc(g->n_mgau, sizeof(int8));

    mg = (ps_mgau_t *)msg;
    mg->vt = &ms_mgau_funcs;
    return mg;
error_out:
    ms_mgau_free(ps_mgau_base(msg));
    return NULL;    
}

void
ms_mgau_free(ps_mgau_t * mg)
{
    ms_mgau_model_t *msg = (ms_mgau_model_t *)mg;
    if (msg == NULL)
        return;

    if (msg->g)
	gauden_free(msg->g);
    if (msg->s)
        senone_free(msg->s);
    if (msg->dist)
        ckd_free_3d((void *) msg->dist);
    if (msg->mgau_active)
        ckd_free(msg->mgau_active);
    
    ckd_free(msg);
}

int
ms_mgau_mllr_transform(ps_mgau_t *s,
		       ps_mllr_t *mllr)
{
    ms_mgau_model_t *msg = (ms_mgau_model_t *)s;
    return gauden_mllr_transform(msg->g, mllr, msg->config);
}

int32
ms_cont_mgau_frame_eval(ps_mgau_t * mg,
			int16 *senscr,
			uint8 *senone_active,
			int32 n_senone_active,
                        mfcc_t ** feat,
			int32 frame,
			int32 compallsen)
{
    ms_mgau_model_t *msg = (ms_mgau_model_t *)mg;
    int32 gid;
    int32 topn;
    int32 best;
    gauden_t *g;
    senone_t *sen;

    topn = ms_mgau_topn(msg);
    g = ms_mgau_gauden(msg);
    sen = ms_mgau_senone(msg);

    if (compallsen) {
	int32 s;

	for (gid = 0; gid < g->n_mgau; gid++)
	    gauden_dist(g, gid, topn, feat, msg->dist[gid]);

	best = (int32) 0x7fffffff;
	for (s = 0; s < sen->n_sen; s++) {
	    senscr[s] = senone_eval(sen, s, msg->dist[sen->mgau[s]], topn);
	    if (best > senscr[s]) {
		best = senscr[s];
	    }
	}

	/* Normalize senone scores */
	for (s = 0; s < sen->n_sen; s++) {
	    int32 bs = senscr[s] - best;
	    if (bs > 32767)
		bs = 32767;
	    if (bs < -32768)
		bs = -32768;
	    senscr[s] = bs;
	}
    }
    else {
	int32 i, n;
	/* Flag all active mixture-gaussian codebooks */
	for (gid = 0; gid < g->n_mgau; gid++)
	    msg->mgau_active[gid] = 0;

	n = 0;
	for (i = 0; i < n_senone_active; i++) {
	    /* senone_active consists of deltas. */
	    int32 s = senone_active[i] + n;
	    msg->mgau_active[sen->mgau[s]] = 1;
	    n = s;
	}

	/* Compute topn gaussian density values (for active codebooks) */
	for (gid = 0; gid < g->n_mgau; gid++) {
	    if (msg->mgau_active[gid])
		gauden_dist(g, gid, topn, feat, msg->dist[gid]);
	}

	best = (int32) 0x7fffffff;
	n = 0;
	for (i = 0; i < n_senone_active; i++) {
	    int32 s = senone_active[i] + n;
	    senscr[s] = senone_eval(sen, s, msg->dist[sen->mgau[s]], topn);
	    if (best > senscr[s]) {
		best = senscr[s];
	    }
	    n = s;
	}

	/* Normalize senone scores */
	n = 0;
	for (i = 0; i < n_senone_active; i++) {
	    int32 s = senone_active[i] + n;
	    int32 bs = senscr[s] - best;
	    if (bs > 32767)
		bs = 32767;
	    if (bs < -32768)
		bs = -32768;
	    senscr[s] = bs;
	    n = s;
	}
    }

    return 0;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * ms_mgau.h -- Essentially a wrapper that wrap up gauden and
 * senone. It supports multi-stream. 
 *
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1997 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * HISTORY
 * $Log$
 * Revision 1.1  2006/04/05  20:27:30  dhdfu
 * A Great Reorganzation of header files and executables
 * 
 * Revision 1.3  2006/02/22 16:57:15  arthchan2003
 * Fixed minor dox-doc issue
 *
 * Revision 1.2  2006/02/22 16:56:01  arthchan2003
 * Merged from SPHINX3_5_2_RCI_IRII_BRANCH: Added ms_mgau.[ch] into the trunk. It is a wrapper of ms_gauden and ms_senone
 *
 * Revision 1.1.2.4  2005/09/25 18:55:19  arthchan2003
 * Added a flag to turn on and off precomputation.
 *
 * Revision 1.1.2.3  2005/08/03 18:53:44  dhdfu
 * Add memory deallocation functions.  Also move all the initialization
 * of ms_mgau_model_t into ms_mgau_init (duh!), which entails removing it
 * from decode_anytopo and friends.
 *
 * Revision 1.1.2.2  2005/08/02 21:05:38  arthchan2003
 * 1, Added dist and mgau_active as intermediate variable for computation. 2, Added ms_cont_mgau_frame_eval, which is a multi stream version of GMM computation mainly s3.0 family of tools. 3, Fixed dox-doc.
 *
 * Revision 1.1.2.1  2005/07/20 19:37:09  arthchan2003
 * Added a multi-stream cont_mgau (ms_mgau) which is a wrapper of both gauden and senone.  Add ms_mgau_init and model_set_mllr.  This allow eliminating 600 lines of code in decode_anytopo/align/allphone.
 *
 *
 *
 */

/** \file ms_mgau.h
 *
 * \brief (Sphinx 3.0 specific) A module that wraps up the code of
 * gauden and senone because they are closely related.  
 *
 * At the time at Sphinx 3.1 to 3.2, Ravi has decided to rewrite only
 * single-stream part of the code into cont_mgau.[ch].  This marks the
 * beginning of historical problem of having two sets of Gaussian
 * distribution computation routine, one for single-stream and one of
 * multi-stream.
 *
 * In Sphinx 3.5, when we figure out that it is possible to allow both
 * 3.0 family of tools and 3.x family of tools to coexist.  This
 * becomes one problem we found that very hard to reconcile.  That is
 * why we currently allow two versions of the code in the code
 * base. This is likely to change in the future.
 */


#ifndef _LIBFBS_MS_CONT_MGAU_H_
#define _LIBFBS_MS_CONT_MGAU_H_

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/feat.h>

/* Local headers. */
#include "acmod.h"
#include "bin_mdef.h"
#include "ms_gauden.h"
#include "ms_senone.h"

/** \struct ms_mgau_t
    \brief Multi-stream mixture gaussian. It is not necessary to be continr
*/

typedef struct {
    ps_mgau_t base;
    gauden_t* g;   /**< The codebook */
    senone_t* s;   /**< The senone */
    int topn;      /**< Top-n gaussian will be computed */

    /**< Intermediate used in computation */
    gauden_dist_t ***dist;  
    uint8 *mgau_active;
    cmd_ln_t *config;
} ms_mgau_model_t;  

#define ms_mgau_gauden(msg) (msg->g)
#define ms_mgau_senone(msg) (msg->s)
#define ms_mgau_topn(msg) (msg->topn)

ps_mgau_t* ms_mgau_init(acmod_t *acmod, logmath_t *lmath, bin_mdef_t *mdef);
void ms_mgau_free(ps_mgau_t *g);
int32 ms_cont_mgau_frame_eval(ps_mgau_t * msg,
                              int16 *senscr,
                              uint8 *senone_active,
                              int32 n_senone_active,
                              mfcc_t ** feat,
                              int32 frame,
                              int32 compallsen);
int32 ms_mgau_mllr_transform(ps_mgau_t *s,
                             ps_mllr_t *mllr);

#endif /* _LIBFBS_MS_CONT_MGAU_H_*/

/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/* System headers. */
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/bio.h>

/* Local headers. */
#include "ms_senone.h"

#define MIXW_PARAM_VERSION	"1.0"
#define SPDEF_PARAM_VERSION	"1.2"

static int32
senone_mgau_map_read(senone_t * s, char const *file_name)
{
    FILE *fp;
    int32 byteswap, chksum_present, n_gauden_present;
    uint32 chksum;
    int32 i;
    char eofchk;
    char **argname, **argval;
    void *ptr;
    float32 v;

    E_INFO("Reading senone gauden-codebook map file: %s\n", file_name);

    if ((fp = fopen(file_name, "rb")) == NULL)
        E_FATAL_SYSTEM("Failed to open map file '%s' for reading", file_name);

    /* Read header, including argument-value info and 32-bit byteorder magic */
    if (bio_readhdr(fp, &argname, &argval, &byteswap) < 0)
        E_FATAL("Failed to read header from file '%s'\n", file_name);

    /* Parse argument-value list */
    chksum_present = 0;
    n_gauden_present = 0;
    for (i = 0; argname[i]; i++) {
        if (strcmp(argname[i], "version") == 0) {
            if (strcmp(argval[i], SPDEF_PARAM_VERSION) != 0) {
                E_WARN("Version mismatch(%s): %s, expecting %s\n",
                       file_name, argval[i], SPDEF_PARAM_VERSION);
            }

            /* HACK!! Convert version# to float32 and take appropriate action */
            if (sscanf(argval[i], "%f", &v) != 1)
                E_FATAL("%s: Bad version no. string: %s\n", file_name,
                        argval[i]);

            n_gauden_present = (v > 1.1) ? 1 : 0;
        }
        else if (strcmp(argname[i], "chksum0") == 0) {
            chksum_present = 1; /* Ignore the associated value */
        }
    }
    bio_hdrarg_free(argname, argval);
    argname = argval = NULL;

    chksum = 0;

    /* Read #gauden (if version matches) */
    if (n_gauden_present) {
        E_INFO("Reading number of codebooks from %s\n", file_name);
        if (bio_fread
            (&(s->n_gauden), sizeof(int32), 1, fp, byteswap, &chksum) != 1)
            E_FATAL("fread(%s) (#gauden) failed\n", file_name);
    }

    /* Read 1d array data */
    if (bio_fread_1d(&ptr, sizeof(uint32), &(s->n_sen), fp,
		     byteswap, &chksum) < 0) {
        E_FATAL("bio_fread_1d(%s) failed\n", file_name);
    }
    s->mgau = ptr;
    E_INFO("Mapping %d senones to %d codebooks\n", s->n_sen, s->n_gauden);

    /* Infer n_gauden if not present in this version */
    if (!n_gauden_present) {
        s->n_gauden = 1;
        for (i = 0; i < s->n_sen; i++)
            if (s->mgau[i] >= s->n_gauden)
                s->n_gauden = s->mgau[i] + 1;
    }

    if (chksum_present)
        bio_verify_chksum(fp, byteswap, chksum);

    if (fread(&eofchk, 1, 1, fp) == 1)
        E_FATAL("More data than expected in %s: %d\n", file_name, eofchk);

    fclose(fp);

    E_INFO("Read %d->%d senone-codebook mappings\n", s->n_sen,
           s->n_gauden);

    return 1;
}


static int32
senone_mixw_read(senone_t * s, char const *file_name, logmath_t *lmath)
{
    char eofchk;
    FILE *fp;
    int32 byteswap, chksum_present;
    uint32 chksum;
    float32 *pdf;
    int32 i, f, c, p, n_err;
    char **argname, **argval;

    E_INFO("Reading senone mixture weights: %s\n", file_name);

    if ((fp = fopen(file_name, "rb")) == NULL)
        E_FATAL_SYSTEM("Failed to open mixture weights file '%s' for reading", file_name);

    /* Read header, including argument-value info and 32-bit byteorder magic */
    if (bio_readhdr(fp, &argname, &argval, &byteswap) < 0)
        E_FATAL("Failed to read header from file '%s'\n", file_name);

    /* Parse argument-value list */
    chksum_present = 0;
    for (i = 0; argname[i]; i++) {
        if (strcmp(argname[i], "version") == 0) {
            if (strcmp(argval[i], MIXW_PARAM_VERSION) != 0)
                E_WARN("Version mismatch(%s): %s, expecting %s\n",
                       file_name, argval[i], MIXW_PARAM_VERSION);
        }
        else if (strcmp(argname[i], "chksum0") == 0) {
            chksum_present = 1; /* Ignore the associated value */
        }
    }
    bio_hdrarg_free(argname, argval);
    argname = argval = NULL;

    chksum = 0;

    /* Read #senones, #features, #codewords, arraysize */
    if ((bio_fread(&(s->n_sen), sizeof(int32), 1, fp, byteswap, &chksum) !=
         1)
        ||
        (bio_fread(&(s->n_feat), sizeof(int32), 1, fp, byteswap, &chksum)
         != 1)
        || (bio_fread(&(s->n_cw), sizeof(int32), 1, fp, byteswap, &chksum)
            != 1)
        || (bio_fread(&i, sizeof(int32), 1, fp, byteswap, &chksum) != 1)) {
        E_FATAL("bio_fread(%s) (arraysize) failed\n", file_name);
    }
    if (i != s->n_sen * s->n_feat * s->n_cw) {
        E_FATAL
            ("%s: #float32s(%d) doesn't match dimensions: %d x %d x %d\n",
             file_name, i, s->n_sen, s->n_feat, s->n_cw);
    }

    /*
     * Compute #LSB bits to be dropped to represent mixwfloor with 8 bits.
     * All PDF values will be truncated (in the LSB positions) by these many bits.
     */
    if ((s->mixwfloor <= 0.0) || (s->mixwfloor >= 1.0))
        E_FATAL("mixwfloor (%e) not in range (0, 1)\n", s->mixwfloor);

    /* Use a fixed shift for compatibility with everything else. */
    E_INFO("Truncating senone logs3(pdf) values by %d bits\n", SENSCR_SHIFT);

    /*
     * Allocate memory for senone PDF data.  Organize normally or transposed depending on
     * s->n_gauden.
     */
    if (s->n_gauden > 1) {
	E_INFO("Not transposing mixture weights in memory\n");
        s->pdf =
            (senprob_t ***) ckd_calloc_3d(s->n_sen, s->n_feat, s->n_cw,
                                          sizeof(senprob_t));
    }
    else {
	E_INFO("Transposing mixture weights in memory\n");
        s->pdf =
            (senprob_t ***) ckd_calloc_3d(s->n_feat, s->n_cw, s->n_sen,
                                          sizeof(senprob_t));
    }

    /* Temporary structure to read in floats */
    pdf = (float32 *) ckd_calloc(s->n_cw, sizeof(float32));

    /* Read senone probs data, normalize, floor, convert to logs3, truncate to 8 bits */
    n_err = 0;
    for (i = 0; i < s->n_sen; i++) {
        for (f = 0; f < s->n_feat; f++) {
            if (bio_fread
                ((void *) pdf, sizeof(float32), s->n_cw, fp, byteswap,
                 &chksum)
                != s->n_cw) {
                E_FATAL("bio_fread(%s) (arraydata) failed\n", file_name);
            }

            /* Normalize and floor */
            if (vector_sum_norm(pdf, s->n_cw) <= 0.0)
                n_err++;
            vector_floor(pdf, s->n_cw, s->mixwfloor);
            vector_sum_norm(pdf, s->n_cw);

            /* Convert to logs3, truncate to 8 bits, and store in s->pdf */
            for (c = 0; c < s->n_cw; c++) {
                p = -(logmath_log(lmath, pdf[c]));
                p += (1 << (SENSCR_SHIFT - 1)) - 1; /* Rounding before truncation */

                if (s->n_gauden > 1)
                    s->pdf[i][f][c] =
                        (p < (255 << SENSCR_SHIFT)) ? (p >> SENSCR_SHIFT) : 255;
                else
                    s->pdf[f][c][i] =
                        (p < (255 << SENSCR_SHIFT)) ? (p >> SENSCR_SHIFT) : 255;
            }
        }
    }
    if (n_err > 0)
        E_WARN("Weight normalization failed for %d mixture weights components\n", n_err);

    ckd_free(pdf);

    if (chksum_present)
        bio_verify_chksum(fp, byteswap, chksum);

    if (fread(&eofchk, 1, 1, fp) == 1)
        E_FATAL("More data than expected in %s\n", file_name);

    fclose(fp);

    E_INFO
        ("Read mixture weights for %d senones: %d features x %d codewords\n",
         s->n_sen, s->n_feat, s->n_cw);

    return 1;
}


senone_t *
senone_init(gauden_t *g, char const *mixwfile, char const *sen2mgau_map_file,
	    float32 mixwfloor, logmath_t *lmath, bin_mdef_t *mdef)
{
    senone_t *s;
    int32 n = 0, i;

    s = (senone_t *) ckd_calloc(1, sizeof(senone_t));
    s->lmath = logmath_init(logmath_get_base(lmath), SENSCR_SHIFT, TRUE);
    s->mixwfloor = mixwfloor;

    s->n_gauden = g->n_mgau;
    if (sen2mgau_map_file) {
	if (!(strcmp(sen2mgau_map_file, ".semi.") == 0
	      || strcmp(sen2mgau_map_file, ".ptm.") == 0
	      || strcmp(sen2mgau_map_file, ".cont.") == 0)) {
	    senone_mgau_map_read(s, sen2mgau_map_file);
	    n = s->n_sen;
	}
    }
    else {
	if (s->n_gauden == 1)
	    sen2mgau_map_file = ".semi.";
	else if (s->n_gauden == bin_mdef_n_ciphone(mdef))
	    sen2mgau_map_file = ".ptm.";
	else
	    sen2mgau_map_file = ".cont.";
    }

    senone_mixw_read(s, mixwfile, lmath);

    if (strcmp(sen2mgau_map_file, ".semi.") == 0) {
        /* All-to-1 senones-codebook mapping */
	E_INFO("Mapping all senones to one codebook\n");
        s->mgau = (uint32 *) ckd_calloc(s->n_sen, sizeof(*s->mgau));
    }
    else if (strcmp(sen2mgau_map_file, ".ptm.") == 0) {
        /* All-to-ciphone-id senones-codebook mapping */
	E_INFO("Mapping senones to context-independent phone codebooks\n");
        s->mgau = (uint32 *) ckd_calloc(s->n_sen, sizeof(*s->mgau));
        for (i = 0; i < s->n_sen; i++)
	    s->mgau[i] = bin_mdef_sen2cimap(mdef, i);
    }
    else if (strcmp(sen2mgau_map_file, ".cont.") == 0
             || strcmp(sen2mgau_map_file, ".s3cont.") == 0) {
        /* 1-to-1 senone-codebook mapping */
	E_INFO("Mapping senones to individual codebooks\n");
        if (s->n_sen <= 1)
            E_FATAL("#senone=%d; must be >1\n", s->n_sen);

        s->mgau = (uint32 *) ckd_calloc(s->n_sen, sizeof(*s->mgau));
        for (i = 0; i < s->n_sen; i++)
            s->mgau[i] = i;
	/* Not sure why this is here, it probably does nothing. */
        s->n_gauden = s->n_sen;
    }
    else {
        if (s->n_sen != n)
            E_FATAL("#senones inconsistent: %d in %s; %d in %s\n",
                    n, sen2mgau_map_file, s->n_sen, mixwfile);
    }

    s->featscr = NULL;
    return s;
}

void
senone_free(senone_t * s)
{
    if (s == NULL)
        return;
    if (s->pdf)
        ckd_free_3d((void *) s->pdf);
    if (s->mgau)
        ckd_free(s->mgau);
    if (s->featscr)
        ckd_free(s->featscr);
    logmath_free(s->lmath);
    ckd_free(s);
}


/*
 * Compute senone score for one senone.
 * NOTE:  Remember that senone PDF tables contain SCALED, NEGATED logs3 values.
 * NOTE:  Remember also that PDF data may be transposed or not depending on s->n_gauden.
 */
int32
senone_eval(senone_t * s, int id, gauden_dist_t ** dist, int32 n_top)
{
    int32 scr;                  /* total senone score */
    int32 fden;                 /* Gaussian density */
    int32 fscr;                 /* senone score for one feature */
    int32 fwscr;                /* senone score for one feature, one codeword */
    int32 f, t;
    gauden_dist_t *fdist;

    assert((id >= 0) && (id < s->n_sen));
    assert((n_top > 0) && (n_top <= s->n_cw));

    scr = 0;

    for (f = 0; f < s->n_feat; f++) {
        int top;
        fdist = dist[f];

        /* Top codeword for feature f */
	top = fden = ((int32)fdist[0].dist + ((1<<SENSCR_SHIFT) - 1)) >> SENSCR_SHIFT;
        fscr = (s->n_gauden > 1)
	    ? (fden + -s->pdf[id][f][fdist[0].id])  /* untransposed */
	    : (fden + -s->pdf[f][fdist[0].id][id]); /* transposed */
        E_DEBUG(1, ("fden[%d][%d] l+= %d + %d = %d\n",
                    id, f, -(fscr - fden), -(fden-top), -(fscr-top)));
        /* Remaining of n_top codewords for feature f */
        for (t = 1; t < n_top; t++) {
	    fden = ((int32)fdist[t].dist + ((1<<SENSCR_SHIFT) - 1)) >> SENSCR_SHIFT;
            fwscr = (s->n_gauden > 1) ?
                (fden + -s->pdf[id][f][fdist[t].id]) :
                (fden + -s->pdf[f][fdist[t].id][id]);
            fscr = logmath_add(s->lmath, fscr, fwscr);
            E_DEBUG(1, ("fden[%d][%d] l+= %d + %d = %d\n",
                        id, f, -(fwscr - fden), -(fden-top), -(fscr-top)));
        }
	/* Senone scores are also scaled, negated logs3 values.  Hence
	 * we have to negate the stuff we calculated above. */
        scr -= fscr;
    }
    /* Downscale scores. */
    scr /= s->aw;

    /* Avoid overflowing int16 */
    if (scr > 32767)
      scr = 32767;
    if (scr < -32768)
      scr = -32768;
    return scr;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * senone.h -- Mixture density weights associated with each tied state.
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1996 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * 
 * $Log$
 * Revision 1.1  2006/04/05  20:27:30  dhdfu
 * A Great Reorganzation of header files and executables
 * 
 * Revision 1.7  2006/02/22 17:27:39  arthchan2003
 * Merged from SPHINX3_5_2_RCI_IRII_BRANCH: 1, NOT doing truncation in the multi-stream GMM computation \n. 2, Added .s3cont. to be the alias of the old multi-stream GMM computation routine \n. 3, Added license \n.  4, Fixed dox-doc. \n
 *
 * Revision 1.6.4.4  2006/01/16 19:47:05  arthchan2003
 * Removed the truncation of senone probability code.
 *
 * Revision 1.6.4.3  2005/08/03 18:53:43  dhdfu
 * Add memory deallocation functions.  Also move all the initialization
 * of ms_mgau_model_t into ms_mgau_init (duh!), which entails removing it
 * from decode_anytopo and friends.
 *
 * Revision 1.6.4.2  2005/07/20 19:39:01  arthchan2003
 * Added licences in ms_* series of code.
 *
 * Revision 1.6.4.1  2005/07/05 05:47:59  arthchan2003
 * Fixed dox-doc. struct level of documentation are included.
 *
 * Revision 1.6  2005/06/21 19:00:19  arthchan2003
 * Add more detail comments  to ms_senone.h
 *
 * Revision 1.5  2005/06/21 18:57:31  arthchan2003
 * 1, Fixed doxygen documentation. 2, Added $ keyword.
 *
 * Revision 1.2  2005/06/13 04:02:56  archan
 * Fixed most doxygen-style documentation under libs3decoder.
 *
 * Revision 1.1.1.1  2005/03/24 15:24:00  archan
 * I found Evandro's suggestion is quite right after yelling at him 2 days later. So I decide to check this in again without any binaries. (I have done make distcheck. ) . Again, this is a candidate for s3.6 and I believe I need to work out 4-5 intermediate steps before I can complete the first prototype.  That's why I keep local copies. 
 *
 * Revision 1.4  2004/12/06 10:52:01  arthchan2003
 * Enable doxygen documentation in libs3decoder
 *
 * Revision 1.3  2004/11/13 21:25:19  arthchan2003
 * commit of 1, absolute CI-GMMS , 2, fast CI senone computation using svq, 3, Decrease the number of static variables, 4, fixing the random generator problem of vector_vqgen, 5, move all unused files to NOTUSED
 *
 * Revision 1.2  2004/08/31 08:43:47  arthchan2003
 * Fixing _cpluscplus directive
 *
 * Revision 1.1  2004/08/09 00:17:11  arthchan2003
 * Incorporating s3.0 align, at this point, there are still some small problems in align but they don't hurt. For example, the score doesn't match with s3.0 and the output will have problem if files are piped to /dev/null/. I think we can go for it.
 *
 * Revision 1.1  2003/02/14 14:40:34  cbq
 * Compiles.  Analysis is probably hosed.
 *
 * Revision 1.1  2000/04/24 09:39:41  lenzo
 * s3 import.
 *
 * 
 * 13-Dec-95	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added senone_eval_all().
 * 
 * 12-Nov-95	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Created.
 */


#ifndef _LIBFBS_SENONE_H_
#define _LIBFBS_SENONE_H_


/* SphinxBase headers. */
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>

/* Local headers. */
#include "ms_gauden.h"
#include "bin_mdef.h"

/** \file ms_senone.h
 *  \brief (Sphinx 3.0 specific) multiple streams senones. used with ms_gauden.h
 *  In Sphinx 3.0 family of tools, ms_senone is used to combine the Gaussian scores.
 *  Its existence is crucial in Sphinx 3.0 because 3.0 supports both SCHMM and CDHMM. 
 *  There are optimization scheme for SCHMM (e.g. compute the top-N Gaussian) that is 
 *  applicable to SCHMM than CDHMM.  This is wrapped in senone_eval_all. 
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

typedef uint8 senprob_t;	/**< Senone logs3-probs, truncated to 8 bits */

/**
 * \struct senone_t
 * \brief 8-bit senone PDF structure. 
 * 
 * 8-bit senone PDF structure.  Senone pdf values are normalized, floored, converted to
 * logs3 domain, and finally truncated to 8 bits precision to conserve memory space.
 */
typedef struct {
    senprob_t ***pdf;		/**< gaussian density mixture weights, organized two possible
                                   ways depending on n_gauden:
                                   if (n_gauden > 1): pdf[sen][feat][codeword].  Not an
                                   efficient representation--memory access-wise--but
                                   evaluating the many codebooks will be more costly.
                                   if (n_gauden == 1): pdf[feat][codeword][sen].  Optimized
                                   for the shared-distribution semi-continuous case. */
    logmath_t *lmath;           /**< log math computation */
    uint32 n_sen;		/**< Number senones in this set */
    uint32 n_feat;		/**< Number feature streams */ 
    uint32 n_cw;		/**< Number codewords per codebook,stream */
    uint32 n_gauden;		/**< Number gaussian density codebooks referred to by senones */
    float32 mixwfloor;		/**< floor applied to each PDF entry */
    uint32 *mgau;		/**< senone-id -> mgau-id mapping for senones in this set */
    int32 *featscr;              /**< The feature score for every senone, will be initialized inside senone_eval_all */
    int32 aw;			/**< Inverse acoustic weight */
} senone_t;


/**
 * Load a set of senones (mixing weights and mixture gaussian codebook mappings) from
 * the given files.  Normalize weights for each codebook, apply the given floor, convert
 * PDF values to logs3 domain and quantize to 8-bits.
 * @return pointer to senone structure created.  Caller MUST NOT change its contents.
 */
senone_t *senone_init (gauden_t *g,             /**< In: codebooks */
                       char const *mixwfile,	/**< In: mixing weights file */
		       char const *mgau_mapfile,/**< In: file specifying mapping from each
						   senone to mixture gaussian codebook.
						   If NULL all senones map to codebook 0 */
		       float32 mixwfloor,	/**< In: Floor value for senone weights */
                       logmath_t *lmath,        /**< In: log math computation */
                       bin_mdef_t *mdef         /**< In: model definition */
    );

/** Release memory allocated by senone_init. */
void senone_free(senone_t *s); /**< In: The senone_t to free */

/**
 * Evaluate the score for the given senone wrt to the given top N gaussian codewords.
 * @return senone score (in logs3 domain).
 */
int32 senone_eval (senone_t *s, int id,		/**< In: senone for which score desired */
		   gauden_dist_t **dist,	/**< In: top N codewords and densities for
						   all features, to be combined into
						   senone score.  IE, dist[f][i] = i-th
						   best <codeword,density> for feaure f */
		   int n_top		/**< In: Length of dist[f], for each f */
    );

#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif

#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ngram_search.c N-Gram based multi-pass search ("FBS")
 */

/* System headers. */
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "pocketsphinx_internal.h"
#include "ps_lattice_internal.h"
#include "ngram_search.h"
#include "ngram_search_fwdtree.h"
#include "ngram_search_fwdflat.h"

static int ngram_search_start(ps_search_t *search);
static int ngram_search_step(ps_search_t *search, int frame_idx);
static int ngram_search_finish(ps_search_t *search);
static int ngram_search_reinit(ps_search_t *search, dict_t *dict, dict2pid_t *d2p);
static char const *ngram_search_hyp(ps_search_t *search, int32 *out_score, int32 *out_is_final);
static int32 ngram_search_prob(ps_search_t *search);
static ps_seg_t *ngram_search_seg_iter(ps_search_t *search, int32 *out_score);

static ps_searchfuncs_t ngram_funcs = {
    /* name: */   "ngram",
    /* start: */  ngram_search_start,
    /* step: */   ngram_search_step,
    /* finish: */ ngram_search_finish,
    /* reinit: */ ngram_search_reinit,
    /* free: */   ngram_search_free,
    /* lattice: */  ngram_search_lattice,
    /* hyp: */      ngram_search_hyp,
    /* prob: */     ngram_search_prob,
    /* seg_iter: */ ngram_search_seg_iter,
};

static void
ngram_search_update_widmap(ngram_search_t *ngs)
{
    const char **words;
    int32 i, n_words;

    /* It's okay to include fillers since they won't be in the LM */
    n_words = ps_search_n_words(ngs);
    words = ckd_calloc(n_words, sizeof(*words));
    /* This will include alternates, again, that's okay since they aren't in the LM */
    for (i = 0; i < n_words; ++i)
        words[i] = (const char *)dict_wordstr(ps_search_dict(ngs), i);
    ngram_model_set_map_words(ngs->lmset, words, n_words);
    ckd_free(words);
}

static void
ngram_search_calc_beams(ngram_search_t *ngs)
{
    cmd_ln_t *config;
    acmod_t *acmod;

    config = ps_search_config(ngs);
    acmod = ps_search_acmod(ngs);

    /* Log beam widths. */
    ngs->beam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-beam"))>>SENSCR_SHIFT;
    ngs->wbeam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-wbeam"))>>SENSCR_SHIFT;
    ngs->pbeam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-pbeam"))>>SENSCR_SHIFT;
    ngs->lpbeam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-lpbeam"))>>SENSCR_SHIFT;
    ngs->lponlybeam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-lponlybeam"))>>SENSCR_SHIFT;
    ngs->fwdflatbeam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-fwdflatbeam"))>>SENSCR_SHIFT;
    ngs->fwdflatwbeam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-fwdflatwbeam"))>>SENSCR_SHIFT;

    /* Absolute pruning parameters. */
    ngs->maxwpf = cmd_ln_int32_r(config, "-maxwpf");
    ngs->maxhmmpf = cmd_ln_int32_r(config, "-maxhmmpf");

    /* Various penalties which may or may not be useful. */
    ngs->wip = logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-wip")) >>SENSCR_SHIFT;
    ngs->nwpen = logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-nwpen")) >>SENSCR_SHIFT;
    ngs->pip = logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-pip")) >>SENSCR_SHIFT;
    ngs->silpen = ngs->pip
        + (logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-silprob"))>>SENSCR_SHIFT);
    ngs->fillpen = ngs->pip
        + (logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-fillprob"))>>SENSCR_SHIFT);

    /* Language weight ratios for fwdflat and bestpath search. */
    ngs->fwdflat_fwdtree_lw_ratio =
        cmd_ln_float32_r(config, "-fwdflatlw")
        / cmd_ln_float32_r(config, "-lw");
    ngs->bestpath_fwdtree_lw_ratio =
        cmd_ln_float32_r(config, "-bestpathlw")
        / cmd_ln_float32_r(config, "-lw");

    /* Acoustic score scale for posterior probabilities. */
    ngs->ascale = 1.0 / cmd_ln_float32_r(config, "-ascale");
}

ps_search_t *
ngram_search_init(cmd_ln_t *config,
		  acmod_t *acmod,
		  dict_t *dict,
                  dict2pid_t *d2p)
{
    ngram_search_t *ngs;
    const char *path;

    ngs = ckd_calloc(1, sizeof(*ngs));
    ps_search_init(&ngs->base, &ngram_funcs, config, acmod, dict, d2p);
    ngs->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
                                   acmod->tmat->tp, NULL, acmod->mdef->sseq);
    if (ngs->hmmctx == NULL) {
        ps_search_free(ps_search_base(ngs));
        return NULL;
    }
    ngs->chan_alloc = listelem_alloc_init(sizeof(chan_t));
    ngs->root_chan_alloc = listelem_alloc_init(sizeof(root_chan_t));
    ngs->latnode_alloc = listelem_alloc_init(sizeof(ps_latnode_t));

    /* Calculate various beam widths and such. */
    ngram_search_calc_beams(ngs);

    /* Allocate a billion different tables for stuff. */
    ngs->word_chan = ckd_calloc(dict_size(dict),
                                sizeof(*ngs->word_chan));
    ngs->word_lat_idx = ckd_calloc(dict_size(dict),
                                   sizeof(*ngs->word_lat_idx));
    ngs->word_active = bitvec_alloc(dict_size(dict));
    ngs->last_ltrans = ckd_calloc(dict_size(dict),
                                  sizeof(*ngs->last_ltrans));

    /* FIXME: All these structures need to be made dynamic with
     * garbage collection. */
    ngs->bp_table_size = cmd_ln_int32_r(config, "-latsize");
    ngs->bp_table = ckd_calloc(ngs->bp_table_size,
                               sizeof(*ngs->bp_table));
    /* FIXME: This thing is frickin' huge. */
    ngs->bscore_stack_size = ngs->bp_table_size * 20;
    ngs->bscore_stack = ckd_calloc(ngs->bscore_stack_size,
                                   sizeof(*ngs->bscore_stack));
    ngs->n_frame_alloc = 256;
    ngs->bp_table_idx = ckd_calloc(ngs->n_frame_alloc + 1,
                                   sizeof(*ngs->bp_table_idx));
    ++ngs->bp_table_idx; /* Make bptableidx[-1] valid */

    /* Allocate active word list array */
    ngs->active_word_list = ckd_calloc_2d(2, dict_size(dict),
                                          sizeof(**ngs->active_word_list));

    /* Load language model(s) */
    if ((path = cmd_ln_str_r(config, "-lmctl"))) {
        ngs->lmset = ngram_model_set_read(config, path, acmod->lmath);
        if (ngs->lmset == NULL) {
            E_ERROR("Failed to read language model control file: %s\n",
                    path);
            goto error_out;
        }
        /* Set the default language model if needed. */
        if ((path = cmd_ln_str_r(config, "-lmname"))) {
            ngram_model_set_select(ngs->lmset, path);
        }
    }
    else if ((path = cmd_ln_str_r(config, "-lm"))) {
        static const char *name = "default";
        ngram_model_t *lm;

        lm = ngram_model_read(config, path, NGRAM_AUTO, acmod->lmath);
        if (lm == NULL) {
            E_ERROR("Failed to read language model file: %s\n", path);
            goto error_out;
        }
        ngs->lmset = ngram_model_set_init(config,
                                          &lm, (char **)&name,
                                          NULL, 1);
        if (ngs->lmset == NULL) {
            E_ERROR("Failed to initialize language model set\n");
            goto error_out;
        }
    }
    if (ngs->lmset != NULL
        && ngram_wid(ngs->lmset, S3_FINISH_WORD) == ngram_unknown_wid(ngs->lmset)) {
        E_ERROR("Language model/set does not contain </s>, recognition will fail\n");
        goto error_out;
    }

    /* Create word mappings. */
    ngram_search_update_widmap(ngs);

    /* Initialize fwdtree, fwdflat, bestpath modules if necessary. */
    if (cmd_ln_boolean_r(config, "-fwdtree")) {
        ngram_fwdtree_init(ngs);
        ngs->fwdtree = TRUE;
        ngs->fwdtree_perf.name = "fwdtree";
        ptmr_init(&ngs->fwdtree_perf);
    }
    if (cmd_ln_boolean_r(config, "-fwdflat")) {
        ngram_fwdflat_init(ngs);
        ngs->fwdflat = TRUE;
        ngs->fwdflat_perf.name = "fwdflat";
        ptmr_init(&ngs->fwdflat_perf);
    }
    if (cmd_ln_boolean_r(config, "-bestpath")) {
        ngs->bestpath = TRUE;
        ngs->bestpath_perf.name = "bestpath";
        ptmr_init(&ngs->bestpath_perf);
    }

    return (ps_search_t *)ngs;

error_out:
    ngram_search_free((ps_search_t *)ngs);
    return NULL;
}

static int
ngram_search_reinit(ps_search_t *search, dict_t *dict, dict2pid_t *d2p)
{
    ngram_search_t *ngs = (ngram_search_t *)search;
    int old_n_words;
    int rv = 0;

    /* Update the number of words. */
    old_n_words = search->n_words;
    if (old_n_words != dict_size(dict)) {
        search->n_words = dict_size(dict);
        /* Reallocate these temporary arrays. */
        ckd_free(ngs->word_lat_idx);
        ckd_free(ngs->word_active);
        ckd_free(ngs->last_ltrans);
        ckd_free_2d(ngs->active_word_list);
        ngs->word_lat_idx = ckd_calloc(search->n_words, sizeof(*ngs->word_lat_idx));
        ngs->word_active = bitvec_alloc(search->n_words);
        ngs->last_ltrans = ckd_calloc(search->n_words, sizeof(*ngs->last_ltrans));
        ngs->active_word_list
            = ckd_calloc_2d(2, search->n_words,
                            sizeof(**ngs->active_word_list));
    }

    /* Free old dict2pid, dict */
    ps_search_base_reinit(search, dict, d2p);
    
    if (ngs->lmset == NULL)
	return 0;

    /* Update beam widths. */
    ngram_search_calc_beams(ngs);

    /* Update word mappings. */
    ngram_search_update_widmap(ngs);

    /* Now rebuild lextrees. */
    if (ngs->fwdtree) {
        if ((rv = ngram_fwdtree_reinit(ngs)) < 0)
            return rv;
    }
    if (ngs->fwdflat) {
        if ((rv = ngram_fwdflat_reinit(ngs)) < 0)
            return rv;
    }

    return rv;
}

void
ngram_search_free(ps_search_t *search)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    ps_search_deinit(search);
    if (ngs->fwdtree)
        ngram_fwdtree_deinit(ngs);
    if (ngs->fwdflat)
        ngram_fwdflat_deinit(ngs);
    if (ngs->bestpath) {
        double n_speech = (double)ngs->n_tot_frame
            / cmd_ln_int32_r(ps_search_config(ngs), "-frate");

        E_INFO("TOTAL bestpath %.2f CPU %.3f xRT\n",
               ngs->bestpath_perf.t_tot_cpu,
               ngs->bestpath_perf.t_tot_cpu / n_speech);
        E_INFO("TOTAL bestpath %.2f wall %.3f xRT\n",
               ngs->bestpath_perf.t_tot_elapsed,
               ngs->bestpath_perf.t_tot_elapsed / n_speech);
    }

    hmm_context_free(ngs->hmmctx);
    listelem_alloc_free(ngs->chan_alloc);
    listelem_alloc_free(ngs->root_chan_alloc);
    listelem_alloc_free(ngs->latnode_alloc);
    ngram_model_free(ngs->lmset);

    ckd_free(ngs->word_chan);
    ckd_free(ngs->word_lat_idx);
    bitvec_free(ngs->word_active);
    ckd_free(ngs->bp_table);
    ckd_free(ngs->bscore_stack);
    if (ngs->bp_table_idx != NULL)
        ckd_free(ngs->bp_table_idx - 1);
    ckd_free_2d(ngs->active_word_list);
    ckd_free(ngs->last_ltrans);
    ckd_free(ngs);
}

int
ngram_search_mark_bptable(ngram_search_t *ngs, int frame_idx)
{
    if (frame_idx >= ngs->n_frame_alloc) {
        ngs->n_frame_alloc *= 2;
        ngs->bp_table_idx = ckd_realloc(ngs->bp_table_idx - 1,
                                        (ngs->n_frame_alloc + 1)
                                        * sizeof(*ngs->bp_table_idx));
        if (ngs->frm_wordlist) {
            ngs->frm_wordlist = ckd_realloc(ngs->frm_wordlist,
                                            ngs->n_frame_alloc
                                            * sizeof(*ngs->frm_wordlist));
        }
        ++ngs->bp_table_idx; /* Make bptableidx[-1] valid */
    }
    ngs->bp_table_idx[frame_idx] = ngs->bpidx;
    return ngs->bpidx;
}

static void
set_real_wid(ngram_search_t *ngs, int32 bp)
{
    bptbl_t *ent, *prev;

    assert(bp != NO_BP);
    ent = ngs->bp_table + bp;
    if (ent->bp == NO_BP)
        prev = NULL;
    else
        prev = ngs->bp_table + ent->bp;

    /* Propagate lm state for fillers, rotate it for words. */
    if (dict_filler_word(ps_search_dict(ngs), ent->wid)) {
        if (prev != NULL) {
            ent->real_wid = prev->real_wid;
            ent->prev_real_wid = prev->prev_real_wid;
        }
        else {
            ent->real_wid = dict_basewid(ps_search_dict(ngs),
                                         ent->wid);
            ent->prev_real_wid = BAD_S3WID;
        }
    }
    else {
        ent->real_wid = dict_basewid(ps_search_dict(ngs), ent->wid);
        if (prev != NULL)
            ent->prev_real_wid = prev->real_wid;
        else
            ent->prev_real_wid = BAD_S3WID;
    }
}

#define NGRAM_HISTORY_LONG_WORD 2000 /* 20s */

void
ngram_search_save_bp(ngram_search_t *ngs, int frame_idx,
                     int32 w, int32 score, int32 path, int32 rc)
{
    int32 bp;

    /* Look for an existing exit for this word in this frame.  The
     * only reason one would exist is from a different right context
     * triphone, but of course that happens quite frequently. */
    bp = ngs->word_lat_idx[w];
    if (bp != NO_BP) {

        if (frame_idx - ngs->bp_table[path].frame > NGRAM_HISTORY_LONG_WORD) {
    	    E_WARN("Word '%s' survived for %d frames, potential overpruning\n", dict_wordstr(ps_search_dict(ngs), w),
	    	    frame_idx - ngs->bp_table[path].frame);
	}

        /* Keep only the best scoring one, we will reconstruct the
         * others from the right context scores - usually the history
         * is not lost. */
        if (ngs->bp_table[bp].score WORSE_THAN score) {
            assert(path != bp); /* Pathological. */
            if (ngs->bp_table[bp].bp != path) {
                int32 bplh[2], newlh[2];
                /* But, sometimes, the history *is* lost.  If we wanted to
                 * do exact language model scoring we'd have to preserve
                 * these alternate histories. */
                E_DEBUG(2,("Updating path history %d => %d frame %d\n",
                           ngs->bp_table[bp].bp, path, frame_idx));
                bplh[0] = ngs->bp_table[bp].bp == -1
                    ? -1 : ngs->bp_table[ngs->bp_table[bp].bp].prev_real_wid;
                bplh[1] = ngs->bp_table[bp].bp == -1
                    ? -1 : ngs->bp_table[ngs->bp_table[bp].bp].real_wid;
                newlh[0] = path == -1
                    ? -1 : ngs->bp_table[path].prev_real_wid;
                newlh[1] = path == -1
                    ? -1 : ngs->bp_table[path].real_wid;
                /* Actually it's worth checking how often the actual
                 * language model state changes. */
                if (bplh[0] != newlh[0] || bplh[1] != newlh[1]) {
                    /* It's fairly rare that the actual language model
                     * state changes, but it does happen some
                     * times. */
                    E_DEBUG(1, ("Updating language model state %s,%s => %s,%s frame %d\n",
                                dict_wordstr(ps_search_dict(ngs), bplh[0]),
                                dict_wordstr(ps_search_dict(ngs), bplh[1]),
                                dict_wordstr(ps_search_dict(ngs), newlh[0]),
                                dict_wordstr(ps_search_dict(ngs), newlh[1]),
                                frame_idx));
                    set_real_wid(ngs, bp);
                }
                ngs->bp_table[bp].bp = path;
            }
            ngs->bp_table[bp].score = score;
        }
        /* But do keep track of scores for all right contexts, since
         * we need them to determine the starting path scores for any
         * successors of this word exit. */
        if (ngs->bp_table[bp].s_idx != -1)
            ngs->bscore_stack[ngs->bp_table[bp].s_idx + rc] = score;
    }
    else {
        int32 i, rcsize;
        bptbl_t *be;

        /* This might happen if recognition fails. */
        if (ngs->bpidx == NO_BP) {
            E_ERROR("No entries in backpointer table!");
            return;
        }

        /* Expand the backpointer tables if necessary. */
        if (ngs->bpidx >= ngs->bp_table_size) {
            ngs->bp_table_size *= 2;
            ngs->bp_table = ckd_realloc(ngs->bp_table,
                                        ngs->bp_table_size
                                        * sizeof(*ngs->bp_table));
            E_INFO("Resized backpointer table to %d entries\n", ngs->bp_table_size);
        }
        if (ngs->bss_head >= ngs->bscore_stack_size
            - bin_mdef_n_ciphone(ps_search_acmod(ngs)->mdef)) {
            ngs->bscore_stack_size *= 2;
            ngs->bscore_stack = ckd_realloc(ngs->bscore_stack,
                                            ngs->bscore_stack_size
                                            * sizeof(*ngs->bscore_stack));
            E_INFO("Resized score stack to %d entries\n", ngs->bscore_stack_size);
        }

        ngs->word_lat_idx[w] = ngs->bpidx;
        be = &(ngs->bp_table[ngs->bpidx]);
        be->wid = w;
        be->frame = frame_idx;
        be->bp = path;
        be->score = score;
        be->s_idx = ngs->bss_head;
        be->valid = TRUE;
        assert(path != ngs->bpidx);

        /* DICT2PID */
        /* Get diphone ID for final phone and number of ssids corresponding to it. */
        be->last_phone = dict_last_phone(ps_search_dict(ngs),w);
        if (dict_is_single_phone(ps_search_dict(ngs), w)) {
            be->last2_phone = -1;
            be->s_idx = -1;
            rcsize = 0;
        }
        else {
            be->last2_phone = dict_second_last_phone(ps_search_dict(ngs),w);
            rcsize = dict2pid_rssid(ps_search_dict2pid(ngs),
                                    be->last_phone, be->last2_phone)->n_ssid;
        }
        /* Allocate some space on the bscore_stack for all of these triphones. */
        for (i = 0; i < rcsize; ++i)
            ngs->bscore_stack[ngs->bss_head + i] = WORST_SCORE;
        if (rcsize)
            ngs->bscore_stack[ngs->bss_head + rc] = score;
        set_real_wid(ngs, ngs->bpidx);

        ngs->bpidx++;
        ngs->bss_head += rcsize;
    }
}

int
ngram_search_find_exit(ngram_search_t *ngs, int frame_idx, int32 *out_best_score, int32 *out_is_final)
{
    /* End of backpointers for this frame. */
    int end_bpidx;
    int best_exit, bp;
    int32 best_score;

    /* No hypothesis means no exit node! */
    if (ngs->n_frame == 0)
        return NO_BP;

    if (frame_idx == -1 || frame_idx >= ngs->n_frame)
        frame_idx = ngs->n_frame - 1;
    end_bpidx = ngs->bp_table_idx[frame_idx];

    best_score = WORST_SCORE;
    best_exit = NO_BP;

    /* Scan back to find a frame with some backpointers in it. */
    while (frame_idx >= 0 && ngs->bp_table_idx[frame_idx] == end_bpidx)
        --frame_idx;
    /* This is NOT an error, it just means there is no hypothesis yet. */
    if (frame_idx < 0)
        return NO_BP;

    /* Now find the entry for </s> OR the best scoring entry. */
    assert(end_bpidx < ngs->bp_table_size);
    for (bp = ngs->bp_table_idx[frame_idx]; bp < end_bpidx; ++bp) {
        if (ngs->bp_table[bp].wid == ps_search_finish_wid(ngs)
            || ngs->bp_table[bp].score BETTER_THAN best_score) {
            best_score = ngs->bp_table[bp].score;
            best_exit = bp;
        }
        if (ngs->bp_table[bp].wid == ps_search_finish_wid(ngs))
            break;
    }

    if (out_best_score) {
	*out_best_score = best_score;
    }
    if (out_is_final) {
	*out_is_final = (ngs->bp_table[bp].wid == ps_search_finish_wid(ngs));
    }
    return best_exit;
}

char const *
ngram_search_bp_hyp(ngram_search_t *ngs, int bpidx)
{
    ps_search_t *base = ps_search_base(ngs);
    char *c;
    size_t len;
    int bp;

    if (bpidx == NO_BP)
        return NULL;

    bp = bpidx;
    len = 0;
    while (bp != NO_BP) {
        bptbl_t *be = &ngs->bp_table[bp];
        bp = be->bp;
        if (dict_real_word(ps_search_dict(ngs), be->wid))
            len += strlen(dict_basestr(ps_search_dict(ngs), be->wid)) + 1;
    }

    ckd_free(base->hyp_str);
    if (len == 0) {
	base->hyp_str = NULL;
	return base->hyp_str;
    }
    base->hyp_str = ckd_calloc(1, len);

    bp = bpidx;
    c = base->hyp_str + len - 1;
    while (bp != NO_BP) {
        bptbl_t *be = &ngs->bp_table[bp];
        size_t len;

        bp = be->bp;
        if (dict_real_word(ps_search_dict(ngs), be->wid)) {
            len = strlen(dict_basestr(ps_search_dict(ngs), be->wid));
            c -= len;
            memcpy(c, dict_basestr(ps_search_dict(ngs), be->wid), len);
            if (c > base->hyp_str) {
                --c;
                *c = ' ';
            }
        }
    }

    return base->hyp_str;
}

void
ngram_search_alloc_all_rc(ngram_search_t *ngs, int32 w)
{
    chan_t *hmm, *thmm;
    xwdssid_t *rssid;
    int32 i, tmatid, ciphone;

    /* DICT2PID */
    /* Get pointer to array of triphones for final diphone. */
    assert(!dict_is_single_phone(ps_search_dict(ngs), w));
    ciphone = dict_last_phone(ps_search_dict(ngs),w);
    rssid = dict2pid_rssid(ps_search_dict2pid(ngs),
                           ciphone,
                           dict_second_last_phone(ps_search_dict(ngs),w));
    tmatid = bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef, ciphone);
    hmm = ngs->word_chan[w];
    if ((hmm == NULL) || (hmm_nonmpx_ssid(&hmm->hmm) != rssid->ssid[0])) {
        hmm = listelem_malloc(ngs->chan_alloc);
        hmm->next = ngs->word_chan[w];
        ngs->word_chan[w] = hmm;

        hmm->info.rc_id = 0;
        hmm->ciphone = ciphone;
        hmm_init(ngs->hmmctx, &hmm->hmm, FALSE, rssid->ssid[0], tmatid);
        E_DEBUG(3,("allocated rc_id 0 ssid %d ciphone %d lc %d word %s\n",
                   rssid->ssid[0], hmm->ciphone,
                   dict_second_last_phone(ps_search_dict(ngs),w),
                   dict_wordstr(ps_search_dict(ngs),w)));
    }
    for (i = 1; i < rssid->n_ssid; ++i) {
        if ((hmm->next == NULL) || (hmm_nonmpx_ssid(&hmm->next->hmm) != rssid->ssid[i])) {
            thmm = listelem_malloc(ngs->chan_alloc);
            thmm->next = hmm->next;
            hmm->next = thmm;
            hmm = thmm;

            hmm->info.rc_id = i;
            hmm->ciphone = ciphone;
            hmm_init(ngs->hmmctx, &hmm->hmm, FALSE, rssid->ssid[i], tmatid);
            E_DEBUG(3,("allocated rc_id %d ssid %d ciphone %d lc %d word %s\n",
                       i, rssid->ssid[i], hmm->ciphone,
                       dict_second_last_phone(ps_search_dict(ngs),w),
                       dict_wordstr(ps_search_dict(ngs),w)));
        }
        else
            hmm = hmm->next;
    }
}

void
ngram_search_free_all_rc(ngram_search_t *ngs, int32 w)
{
    chan_t *hmm, *thmm;

    for (hmm = ngs->word_chan[w]; hmm; hmm = thmm) {
        thmm = hmm->next;
        hmm_deinit(&hmm->hmm);
        listelem_free(ngs->chan_alloc, hmm);
    }
    ngs->word_chan[w] = NULL;
}

int32
ngram_search_exit_score(ngram_search_t *ngs, bptbl_t *pbe, int rcphone)
{
    /* DICT2PID */
    /* Get the mapping from right context phone ID to index in the
     * right context table and the bscore_stack. */
    if (pbe->last2_phone == -1) {
        /* No right context for single phone predecessor words. */
        return pbe->score;
    }
    else {
        xwdssid_t *rssid;
        /* Find the index for the last diphone of the previous word +
         * the first phone of the current word. */
        rssid = dict2pid_rssid(ps_search_dict2pid(ngs),
                               pbe->last_phone, pbe->last2_phone);
        /* This may be WORST_SCORE, which means that there was no exit
         * with rcphone as right context. */
        return ngs->bscore_stack[pbe->s_idx + rssid->cimap[rcphone]];
    }
}

/*
 * Compute acoustic and LM scores for a BPTable entry (segment).
 */
void
ngram_compute_seg_score(ngram_search_t *ngs, bptbl_t *be, float32 lwf,
                        int32 *out_ascr, int32 *out_lscr)
{
    bptbl_t *pbe;
    int32 start_score;

    /* Start of utterance. */
    if (be->bp == NO_BP) {
        *out_ascr = be->score;
        *out_lscr = 0;
        return;
    }

    /* Otherwise, calculate lscr and ascr. */
    pbe = ngs->bp_table + be->bp;
    start_score = ngram_search_exit_score(ngs, pbe,
                                 dict_first_phone(ps_search_dict(ngs),be->wid));
    assert(start_score BETTER_THAN WORST_SCORE);

    /* FIXME: These result in positive acoustic scores when filler
       words have non-filler pronunciations.  That whole business
       is still pretty much broken but at least it doesn't
       segfault. */
    if (be->wid == ps_search_silence_wid(ngs)) {
        *out_lscr = ngs->silpen;
    }
    else if (dict_filler_word(ps_search_dict(ngs), be->wid)) {
        *out_lscr = ngs->fillpen;
    }
    else {
        int32 n_used;
        *out_lscr = ngram_tg_score(ngs->lmset,
                                   be->real_wid,
                                   pbe->real_wid,
                                   pbe->prev_real_wid,
                                   &n_used)>>SENSCR_SHIFT;
        *out_lscr = *out_lscr * lwf;
    }
    *out_ascr = be->score - start_score - *out_lscr;
}

static int
ngram_search_start(ps_search_t *search)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    ngs->done = FALSE;
    ngram_model_flush(ngs->lmset);
    if (ngs->fwdtree)
        ngram_fwdtree_start(ngs);
    else if (ngs->fwdflat)
        ngram_fwdflat_start(ngs);
    else
        return -1;
    return 0;
}

static int
ngram_search_step(ps_search_t *search, int frame_idx)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    if (ngs->fwdtree)
        return ngram_fwdtree_search(ngs, frame_idx);
    else if (ngs->fwdflat)
        return ngram_fwdflat_search(ngs, frame_idx);
    else
        return -1;
}

void
dump_bptable(ngram_search_t *ngs)
{
    int i;
    E_INFO("Backpointer table (%d entries):\n", ngs->bpidx);
    for (i = 0; i < ngs->bpidx; ++i) {
        bptbl_t *bpe = ngs->bp_table + i;
        int j, rcsize;

        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d real_wid %-5d prev_real_wid %-5d",
                    i, dict_wordstr(ps_search_dict(ngs), bpe->wid),
                    (bpe->bp == -1
                     ? 0 : ngs->bp_table[bpe->bp].frame + 1),
                    bpe->frame, bpe->score, bpe->bp,
                    bpe->real_wid, bpe->prev_real_wid);

        if (bpe->last2_phone == -1)
            rcsize = 0;
        else
            rcsize = dict2pid_rssid(ps_search_dict2pid(ngs),
                                    bpe->last_phone, bpe->last2_phone)->n_ssid;
        if (rcsize) {
            E_INFOCONT("\tbss");
            for (j = 0; j < rcsize; ++j)
                if (ngs->bscore_stack[bpe->s_idx + j] != WORST_SCORE)
                    E_INFOCONT(" %d", bpe->score - ngs->bscore_stack[bpe->s_idx + j]);
        }
        E_INFOCONT("\n");
    }
}

static int
ngram_search_finish(ps_search_t *search)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    ngs->n_tot_frame += ngs->n_frame;
    if (ngs->fwdtree) {
        ngram_fwdtree_finish(ngs);
        /* dump_bptable(ngs); */

        /* Now do fwdflat search in its entirety, if requested. */
        if (ngs->fwdflat) {
            int i;
            /* Rewind the acoustic model. */
            if (acmod_rewind(ps_search_acmod(ngs)) < 0)
                return -1;
            /* Now redo search. */
            ngram_fwdflat_start(ngs);
            i = 0;
            while (ps_search_acmod(ngs)->n_feat_frame > 0) {
                int nfr;
                if ((nfr = ngram_fwdflat_search(ngs, i)) < 0)
                    return nfr;
                acmod_advance(ps_search_acmod(ngs));
                ++i;
            }
            ngram_fwdflat_finish(ngs);
            /* And now, we should have a result... */
            /* dump_bptable(ngs); */
        }
    }
    else if (ngs->fwdflat) {
        ngram_fwdflat_finish(ngs);
    }

    /* Mark the current utterance as done. */
    ngs->done = TRUE;
    return 0;
}

static ps_latlink_t *
ngram_search_bestpath(ps_search_t *search, int32 *out_score, int backward)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    if (search->last_link == NULL) {
        search->last_link = ps_lattice_bestpath(search->dag, ngs->lmset,
                                                ngs->bestpath_fwdtree_lw_ratio,
                                                ngs->ascale);
        if (search->last_link == NULL)
            return NULL;
        /* Also calculate betas so we can fill in the posterior
         * probability field in the segmentation. */
        if (search->post == 0)
            search->post = ps_lattice_posterior(search->dag, ngs->lmset,
                                                ngs->ascale);
    }
    if (out_score)
        *out_score = search->last_link->path_scr + search->dag->final_node_ascr;
    return search->last_link;
}

static char const *
ngram_search_hyp(ps_search_t *search, int32 *out_score, int32 *out_is_final)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    /* Only do bestpath search if the utterance is complete. */
    if (ngs->bestpath && ngs->done) {
        ps_lattice_t *dag;
        ps_latlink_t *link;
        char const *hyp;
        double n_speech;

        ptmr_reset(&ngs->bestpath_perf);
        ptmr_start(&ngs->bestpath_perf);
        if ((dag = ngram_search_lattice(search)) == NULL)
            return NULL;
        if ((link = ngram_search_bestpath(search, out_score, FALSE)) == NULL)
            return NULL;
        hyp = ps_lattice_hyp(dag, link);
        ptmr_stop(&ngs->bestpath_perf);
        n_speech = (double)dag->n_frames
            / cmd_ln_int32_r(ps_search_config(ngs), "-frate");
        E_INFO("bestpath %.2f CPU %.3f xRT\n",
               ngs->bestpath_perf.t_cpu,
               ngs->bestpath_perf.t_cpu / n_speech);
        E_INFO("bestpath %.2f wall %.3f xRT\n",
               ngs->bestpath_perf.t_elapsed,
               ngs->bestpath_perf.t_elapsed / n_speech);
        return hyp;
    }
    else {
        int32 bpidx;

        /* fwdtree and fwdflat use same backpointer table. */
        bpidx = ngram_search_find_exit(ngs, -1, out_score, out_is_final);
        if (bpidx != NO_BP)
            return ngram_search_bp_hyp(ngs, bpidx);
    }

    return NULL;
}

static void
ngram_search_bp2itor(ps_seg_t *seg, int bp)
{
    ngram_search_t *ngs = (ngram_search_t *)seg->search;
    bptbl_t *be, *pbe;

    be = &ngs->bp_table[bp];
    pbe = be->bp == -1 ? NULL : &ngs->bp_table[be->bp];
    seg->word = dict_wordstr(ps_search_dict(ngs), be->wid);
    seg->ef = be->frame;
    seg->sf = pbe ? pbe->frame + 1 : 0;
    seg->prob = 0; /* Bogus value... */
    /* Compute acoustic and LM scores for this segment. */
    if (pbe == NULL) {
        seg->ascr = be->score;
        seg->lscr = 0;
        seg->lback = 0;
    }
    else {
        int32 start_score;

        /* Find ending path score of previous word. */
        start_score = ngram_search_exit_score(ngs, pbe,
                                     dict_first_phone(ps_search_dict(ngs), be->wid));
        assert(start_score BETTER_THAN WORST_SCORE);
        if (be->wid == ps_search_silence_wid(ngs)) {
            seg->lscr = ngs->silpen;
        }
        else if (dict_filler_word(ps_search_dict(ngs), be->wid)) {
            seg->lscr = ngs->fillpen;
        }
        else {
            seg->lscr = ngram_tg_score(ngs->lmset,
                                       be->real_wid,
                                       pbe->real_wid,
                                       pbe->prev_real_wid,
                                       &seg->lback)>>SENSCR_SHIFT;
            seg->lscr = (int32)(seg->lscr * seg->lwf);
        }
        seg->ascr = be->score - start_score - seg->lscr;
    }
}

static void
ngram_bp_seg_free(ps_seg_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;
    
    ckd_free(itor->bpidx);
    ckd_free(itor);
}

static ps_seg_t *
ngram_bp_seg_next(ps_seg_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;

    if (++itor->cur == itor->n_bpidx) {
        ngram_bp_seg_free(seg);
        return NULL;
    }

    ngram_search_bp2itor(seg, itor->bpidx[itor->cur]);
    return seg;
}

static ps_segfuncs_t ngram_bp_segfuncs = {
    /* seg_next */ ngram_bp_seg_next,
    /* seg_free */ ngram_bp_seg_free
};

static ps_seg_t *
ngram_search_bp_iter(ngram_search_t *ngs, int bpidx, float32 lwf)
{
    bptbl_seg_t *itor;
    int bp, cur;

    /* Calling this an "iterator" is a bit of a misnomer since we have
     * to get the entire backtrace in order to produce it.  On the
     * other hand, all we actually need is the bptbl IDs, and we can
     * allocate a fixed-size array of them. */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &ngram_bp_segfuncs;
    itor->base.search = ps_search_base(ngs);
    itor->base.lwf = lwf;
    itor->n_bpidx = 0;
    bp = bpidx;
    while (bp != NO_BP) {
        bptbl_t *be = &ngs->bp_table[bp];
        bp = be->bp;
        ++itor->n_bpidx;
    }
    if (itor->n_bpidx == 0) {
        ckd_free(itor);
        return NULL;
    }
    itor->bpidx = ckd_calloc(itor->n_bpidx, sizeof(*itor->bpidx));
    cur = itor->n_bpidx - 1;
    bp = bpidx;
    while (bp != NO_BP) {
        bptbl_t *be = &ngs->bp_table[bp];
        itor->bpidx[cur] = bp;
        bp = be->bp;
        --cur;
    }

    /* Fill in relevant fields for first element. */
    ngram_search_bp2itor((ps_seg_t *)itor, itor->bpidx[0]);

    return (ps_seg_t *)itor;
}

static ps_seg_t *
ngram_search_seg_iter(ps_search_t *search, int32 *out_score)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    /* Only do bestpath search if the utterance is done. */
    if (ngs->bestpath && ngs->done) {
        ps_lattice_t *dag;
        ps_latlink_t *link;
        double n_speech;
        ps_seg_t *itor;

        ptmr_reset(&ngs->bestpath_perf);
        ptmr_start(&ngs->bestpath_perf);
        if ((dag = ngram_search_lattice(search)) == NULL)
            return NULL;
        if ((link = ngram_search_bestpath(search, out_score, TRUE)) == NULL)
            return NULL;
        itor = ps_lattice_seg_iter(dag, link,
                                   ngs->bestpath_fwdtree_lw_ratio);
        ptmr_stop(&ngs->bestpath_perf);
        n_speech = (double)dag->n_frames
            / cmd_ln_int32_r(ps_search_config(ngs), "-frate");
        E_INFO("bestpath %.2f CPU %.3f xRT\n",
               ngs->bestpath_perf.t_cpu,
               ngs->bestpath_perf.t_cpu / n_speech);
        E_INFO("bestpath %.2f wall %.3f xRT\n",
               ngs->bestpath_perf.t_elapsed,
               ngs->bestpath_perf.t_elapsed / n_speech);
        return itor;
    }
    else {
        int32 bpidx;

        /* fwdtree and fwdflat use same backpointer table. */
        bpidx = ngram_search_find_exit(ngs, -1, out_score, NULL);
        return ngram_search_bp_iter(ngs, bpidx,
                                    /* but different language weights... */
                                    (ngs->done && ngs->fwdflat)
                                    ? ngs->fwdflat_fwdtree_lw_ratio : 1.0);
    }

    return NULL;
}

static int32
ngram_search_prob(ps_search_t *search)
{
    ngram_search_t *ngs = (ngram_search_t *)search;

    /* Only do bestpath search if the utterance is done. */
    if (ngs->bestpath && ngs->done) {
        ps_lattice_t *dag;
        ps_latlink_t *link;

        if ((dag = ngram_search_lattice(search)) == NULL)
            return 0;
        if ((link = ngram_search_bestpath(search, NULL, TRUE)) == NULL)
            return 0;
        return search->post;
    }
    else {
        /* FIXME: Give some kind of good estimate here, eventually. */
        return 0;
    }
}

static void
create_dag_nodes(ngram_search_t *ngs, ps_lattice_t *dag)
{
    bptbl_t *bp_ptr;
    int32 i;

    for (i = 0, bp_ptr = ngs->bp_table; i < ngs->bpidx; ++i, ++bp_ptr) {
        int32 sf, ef, wid;
        ps_latnode_t *node;

        /* Skip invalid backpointers (these result from -maxwpf pruning) */
        if (!bp_ptr->valid)
            continue;

        sf = (bp_ptr->bp < 0) ? 0 : ngs->bp_table[bp_ptr->bp].frame + 1;
        ef = bp_ptr->frame;
        wid = bp_ptr->wid;

        assert(ef < dag->n_frames);
        /* Skip non-final </s> entries. */
        if ((wid == ps_search_finish_wid(ngs)) && (ef < dag->n_frames - 1))
            continue;

        /* Skip if word not in LM */
        if ((!dict_filler_word(ps_search_dict(ngs), wid))
            && (!ngram_model_set_known_wid(ngs->lmset,
                                           dict_basewid(ps_search_dict(ngs), wid))))
            continue;

        /* See if bptbl entry <wid,sf> already in lattice */
        for (node = dag->nodes; node; node = node->next) {
            if ((node->wid == wid) && (node->sf == sf))
                break;
        }

        /* For the moment, store bptbl indices in node.{fef,lef} */
        if (node)
            node->lef = i;
        else {
            /* New node; link to head of list */
            node = listelem_malloc(dag->latnode_alloc);
            node->wid = wid;
            node->sf = sf; /* This is a frame index. */
            node->fef = node->lef = i; /* These are backpointer indices (argh) */
            node->reachable = FALSE;
            node->entries = NULL;
            node->exits = NULL;

            /* NOTE: This creates the list of nodes in reverse
             * topological order, i.e. a node always precedes its
             * antecedents in this list. */
            node->next = dag->nodes;
            dag->nodes = node;
            ++dag->n_nodes;
        }
    }
}

static ps_latnode_t * 
find_start_node_two(ngram_search_t *ngs, ps_lattice_t *dag)
{
    ps_latnode_t *node;

     /*Find start node <s>.0*/ 
    for (node = dag->nodes; node; node = node->next) {
        if ((node->wid == ps_search_start_wid(ngs)) && (node->sf == 0))
            break;
    }
    if (!node) {
        /* This is probably impossible. */
        E_ERROR("Couldn't find <s> in first frame\n");
        return NULL;
    }
    return node;
}

static ps_latnode_t * find_end_node_two(ngram_search_t *ngs, ps_lattice_t *dag, float32 lwf)
{
    ps_latnode_t *node;
    int32 ef, bestbp, bp, bestscore;

    /* Find final node </s>.last_frame; nothing can follow this node */
    for (node = dag->nodes; node; node = node->next) {
        int32 lef = ngs->bp_table[node->lef].frame;
        if ((node->wid == ps_search_finish_wid(ngs))
            && (lef == dag->n_frames - 1))
            break;
    }
    if (node != NULL)
        return node;

    // It is quite likely that no </s> exited in the last frame.  So,
    // find the node corresponding to the best exit. 
    // Find the last frame containing a word exit. 
    for (ef = dag->n_frames - 1;
         ef >= 0 && ngs->bp_table_idx[ef] == ngs->bpidx;
         --ef);
    if (ef < 0) {
        E_ERROR("Empty backpointer table: can not build DAG.\n");
        return NULL;
    }

     //Find best word exit in that frame. 
    bestscore = WORST_SCORE;
    bestbp = NO_BP;
    for (bp = ngs->bp_table_idx[ef]; bp < ngs->bp_table_idx[ef + 1]; ++bp) {
        int32 n_used, l_scr, wid, prev_wid;
        wid = ngs->bp_table[bp].real_wid;
        prev_wid = ngs->bp_table[bp].prev_real_wid;
        // Always prefer </s>, of which there will only be one per frame. 
        if (wid == ps_search_finish_wid(ngs)) {
            bestbp = bp;
            break;
        }
        l_scr = ngram_tg_score(ngs->lmset, ps_search_finish_wid(ngs),
                               wid, prev_wid, &n_used) >>SENSCR_SHIFT;
        l_scr = l_scr * lwf;
        if (ngs->bp_table[bp].score + l_scr BETTER_THAN bestscore) {
            bestscore = ngs->bp_table[bp].score + l_scr;
            bestbp = bp;
        }
    }
    if (bestbp == NO_BP) {
        E_ERROR("No word exits found in last frame (%d), assuming no recognition\n", ef);
        return NULL;
    }
    E_INFO("</s> not found in last frame, using %s.%d instead\n",
           dict_basestr(ps_search_dict(ngs), ngs->bp_table[bestbp].wid), ef);

    //Now find the node that corresponds to it. 
    for (node = dag->nodes; node; node = node->next) {
        if (node->lef == bestbp)
            return node;
    }

    /* FIXME: This seems to happen a lot! */
    E_ERROR("Failed to find DAG node corresponding to %s\n",
           dict_basestr(ps_search_dict(ngs), ngs->bp_table[bestbp].wid));
    return NULL;
}

/*
 * Build lattice from bptable.
 */
ps_lattice_t *
ngram_search_lattice(ps_search_t *search)
{
    int32 i, score, ascr, lscr;
    ps_latnode_t *node, *from, *to;
    ngram_search_t *ngs;
    ps_lattice_t *dag;
    int min_endfr, nlink;
    float lwf;

    ngs = (ngram_search_t *)search;
    min_endfr = cmd_ln_int32_r(ps_search_config(search), "-min_endfr");

    /* If the best score is WORST_SCORE or worse, there is no way to
     * make a lattice. */
    if (ngs->best_score == WORST_SCORE || ngs->best_score WORSE_THAN WORST_SCORE)
        return NULL;

    /* Check to see if a lattice has previously been created over the
     * same number of frames, and reuse it if so. */
    if (search->dag && search->dag->n_frames == ngs->n_frame)
        return search->dag;

    /* Nope, create a new one. */
    ps_lattice_free(search->dag);
    search->dag = NULL;
    dag = ps_lattice_init_search(search, ngs->n_frame);
    /* Compute these such that they agree with the fwdtree language weight. */
    lwf = ngs->fwdflat ? ngs->fwdflat_fwdtree_lw_ratio : 1.0;
    create_dag_nodes(ngs, dag);
    if ((dag->start = find_start_node_two(ngs, dag)) == NULL)
        goto error_out;
    if ((dag->end = find_end_node_two(ngs, dag, ngs->bestpath_fwdtree_lw_ratio)) == NULL)
        goto error_out;
    E_INFO("lattice start node %s.%d end node %s.%d\n",
           dict_wordstr(search->dict, dag->start->wid), dag->start->sf,
           dict_wordstr(search->dict, dag->end->wid), dag->end->sf);

    ngram_compute_seg_score(ngs, ngs->bp_table + dag->end->lef, lwf,
                            &dag->final_node_ascr, &lscr);

    /*
     * At this point, dag->nodes is ordered such that nodes earlier in
     * the list can follow (in time) those later in the list, but not
     * vice versa (see above - also note that adjacency is purely
     * determined by time which is why we can make this claim).  Now
     * create precedence links and simultanesously mark all nodes that
     * can reach dag->end.  (All nodes are reached from dag->start
     * simply by definition - they were created that way).
     *
     * Note that this also means that any nodes before dag->end in the
     * list can be discarded, meaning that dag->end will always be
     * equal to dag->nodes (FIXME: except when loading from a file but
     * we can fix that...)
     */
    i = 0;
    while (dag->nodes && dag->nodes != dag->end) {
        ps_latnode_t *next = dag->nodes->next;
        listelem_free(dag->latnode_alloc, dag->nodes);
        dag->nodes = next;
        ++i;
    }
    E_INFO("Eliminated %d nodes before end node\n", i);
    dag->end->reachable = TRUE;
    nlink = 0;
    for (to = dag->end; to; to = to->next) {
        int fef, lef;

        /* Skip if not reachable; it will never be reachable from dag->end */
        if (!to->reachable)
            continue;

        /* Prune nodes with too few endpoints - heuristic
           borrowed from Sphinx3 */
        fef = ngs->bp_table[to->fef].frame;
        lef = ngs->bp_table[to->lef].frame;
        if (to != dag->end && lef - fef < min_endfr) {
            to->reachable = FALSE;
            continue;
        }

        /* Find predecessors of to : from->fef+1 <= to->sf <= from->lef+1 */
        for (from = to->next; from; from = from->next) {
            bptbl_t *from_bpe;

            fef = ngs->bp_table[from->fef].frame;
            lef = ngs->bp_table[from->lef].frame;

            if ((to->sf <= fef) || (to->sf > lef + 1))
                continue;
            if (lef - fef < min_endfr) {
                assert(!from->reachable);
                continue;
            }

            /* Find bptable entry for "from" that exactly precedes "to" */
            i = from->fef;
            from_bpe = ngs->bp_table + i;
            for (; i <= from->lef; i++, from_bpe++) {
                if (from_bpe->wid != from->wid)
                    continue;
                if (from_bpe->frame >= to->sf - 1)
                    break;
            }

            if ((i > from->lef) || (from_bpe->frame != to->sf - 1))
                continue;

            /* Find acoustic score from.sf->to.sf-1 with right context = to */
            /* This gives us from_bpe's best acoustic score. */
            ngram_compute_seg_score(ngs, from_bpe, lwf,
                                    &ascr, &lscr);
            /* Now find the exact path score for from->to, including
             * the appropriate final triphone.  In fact this might not
             * exist. */
            score = ngram_search_exit_score(ngs, from_bpe,
                                            dict_first_phone(ps_search_dict(ngs), to->wid));
            /* Does not exist.  Can't create a link here. */
            if (score == WORST_SCORE)
                continue;
            /* Adjust the arc score to match the correct triphone. */
            else
                score = ascr + (score - from_bpe->score);
            if (score BETTER_THAN 0) {
                /* Scores must be negative, or Bad Things will happen.
                   In general, they are, except in corner cases
                   involving filler words.  We don't want to throw any
                   links away so we'll keep these, but with some
                   arbitrarily improbable but recognizable score. */
                ps_lattice_link(dag, from, to, -424242, from_bpe->frame);
                ++nlink;
                from->reachable = TRUE;
            }
            else if (score BETTER_THAN WORST_SCORE) {
                ps_lattice_link(dag, from, to, score, from_bpe->frame);
                ++nlink;
                from->reachable = TRUE;
            }
        }
    }

    /* There must be at least one path between dag->start and dag->end */
    if (!dag->start->reachable) {
        E_ERROR("End node of lattice isolated; unreachable\n");
        goto error_out;
    }

    for (node = dag->nodes; node; node = node->next) {
        /* Change node->{fef,lef} from bptbl indices to frames. */
        node->fef = ngs->bp_table[node->fef].frame;
        node->lef = ngs->bp_table[node->lef].frame;
        /* Find base wid for nodes. */
        node->basewid = dict_basewid(search->dict, node->wid);
    }

    /* Link nodes with alternate pronunciations at the same timepoint. */
    for (node = dag->nodes; node; node = node->next) {
        ps_latnode_t *alt;
        /* Scan forward to find the next alternate, then stop. */
        for (alt = node->next; alt && alt->sf == node->sf; alt = alt->next) {
            if (alt->basewid == node->basewid) {
                alt->alt = node->alt;
                node->alt = alt;
                break;
            }
        }
    }
    E_INFO("Lattice has %d nodes, %d links\n", dag->n_nodes, nlink);

    /* Minor hack: If the final node is a filler word and not </s>,
     * then set its base word ID to </s>, so that the language model
     * scores won't be screwed up. */
    if (dict_filler_word(ps_search_dict(ngs), dag->end->wid))
        dag->end->basewid = ps_search_finish_wid(ngs);

    /* Free nodes unreachable from dag->end and their links */
    ps_lattice_delete_unreachable(dag);

    /* Build links around silence and filler words, since they do not
     * exist in the language model. */
    ps_lattice_bypass_fillers(dag, ngs->silpen, ngs->fillpen);

    search->dag = dag;
    return dag;

error_out:
    ps_lattice_free(dag);
    return NULL;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ngram_search_fwdflat.c Flat lexicon search.
 */

/* System headers. */
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "ngram_search.h"
#include "ps_lattice_internal.h"

/* Turn this on to dump channels for debugging */
#define __CHAN_DUMP__		0
#if __CHAN_DUMP__
#define chan_v_eval(chan) hmm_dump_vit_eval(&(chan)->hmm, stderr)
#else
#define chan_v_eval(chan) hmm_vit_eval(&(chan)->hmm)
#endif

static void
ngram_fwdflat_expand_all(ngram_search_t *ngs)
{
    int n_words, i;

    /* For all "real words" (not fillers or <s>/</s>) in the dictionary,
     *
     * 1) Add the ones which are in the LM to the fwdflat wordlist
     * 2) And to the expansion list (since we are expanding all)
     */
    ngs->n_expand_words = 0;
    n_words = ps_search_n_words(ngs);
    bitvec_clear_all(ngs->expand_word_flag, ps_search_n_words(ngs));
    for (i = 0; i < n_words; ++i) {
        if (!ngram_model_set_known_wid(ngs->lmset,
                                       dict_basewid(ps_search_dict(ngs),i)))
            continue;
        ngs->fwdflat_wordlist[ngs->n_expand_words] = i;
        ngs->expand_word_list[ngs->n_expand_words] = i;
        bitvec_set(ngs->expand_word_flag, i);
        ngs->n_expand_words++;
    }
    E_INFO("Utterance vocabulary contains %d words\n", ngs->n_expand_words);
    ngs->expand_word_list[ngs->n_expand_words] = -1;
    ngs->fwdflat_wordlist[ngs->n_expand_words] = -1;
}

static void
ngram_fwdflat_allocate_1ph(ngram_search_t *ngs)
{
    dict_t *dict = ps_search_dict(ngs);
    int n_words = ps_search_n_words(ngs);
    int i, w;

    /* Allocate single-phone words, since they won't have
     * been allocated for us by fwdtree initialization. */
    ngs->n_1ph_words = 0;
    for (w = 0; w < n_words; w++) {
        if (dict_is_single_phone(dict, w))
            ++ngs->n_1ph_words;
    }
    ngs->single_phone_wid = ckd_calloc(ngs->n_1ph_words,
                                       sizeof(*ngs->single_phone_wid));
    ngs->rhmm_1ph = ckd_calloc(ngs->n_1ph_words, sizeof(*ngs->rhmm_1ph));
    i = 0;
    for (w = 0; w < n_words; w++) {
        if (!dict_is_single_phone(dict, w))
            continue;

        /* DICT2PID location */
        ngs->rhmm_1ph[i].ciphone = dict_first_phone(dict, w);
        ngs->rhmm_1ph[i].ci2phone = bin_mdef_silphone(ps_search_acmod(ngs)->mdef);
        hmm_init(ngs->hmmctx, &ngs->rhmm_1ph[i].hmm, TRUE,
                 /* ssid */ bin_mdef_pid2ssid(ps_search_acmod(ngs)->mdef,
                                              ngs->rhmm_1ph[i].ciphone),
                 /* tmatid */ bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef,
    						  ngs->rhmm_1ph[i].ciphone));
        ngs->rhmm_1ph[i].next = NULL;
        ngs->word_chan[w] = (chan_t *) &(ngs->rhmm_1ph[i]);
        ngs->single_phone_wid[i] = w;
        i++;
    }
}

static void
ngram_fwdflat_free_1ph(ngram_search_t *ngs)
{
    int i, w;
    int n_words = ps_search_n_words(ngs);

    for (i = w = 0; w < n_words; ++w) {
        if (!dict_is_single_phone(ps_search_dict(ngs), w))
            continue;
        hmm_deinit(&ngs->rhmm_1ph[i].hmm);
        ++i;
    }
    ckd_free(ngs->rhmm_1ph);
    ngs->rhmm_1ph = NULL;
    ckd_free(ngs->single_phone_wid);
}

void
ngram_fwdflat_init(ngram_search_t *ngs)
{
    int n_words;

    n_words = ps_search_n_words(ngs);
    ngs->fwdflat_wordlist = ckd_calloc(n_words + 1, sizeof(*ngs->fwdflat_wordlist));
    ngs->expand_word_flag = bitvec_alloc(n_words);
    ngs->expand_word_list = ckd_calloc(n_words + 1, sizeof(*ngs->expand_word_list));
    ngs->frm_wordlist = ckd_calloc(ngs->n_frame_alloc, sizeof(*ngs->frm_wordlist));
    ngs->min_ef_width = cmd_ln_int32_r(ps_search_config(ngs), "-fwdflatefwid");
    ngs->max_sf_win = cmd_ln_int32_r(ps_search_config(ngs), "-fwdflatsfwin");
    E_INFO("fwdflat: min_ef_width = %d, max_sf_win = %d\n",
           ngs->min_ef_width, ngs->max_sf_win);

    /* No tree-search; pre-build the expansion list, including all LM words. */
    if (!ngs->fwdtree) {
        /* Build full expansion list from LM words. */
        ngram_fwdflat_expand_all(ngs);
        /* Allocate single phone words. */
        ngram_fwdflat_allocate_1ph(ngs);
    }
}

void
ngram_fwdflat_deinit(ngram_search_t *ngs)
{
    double n_speech = (double)ngs->n_tot_frame
            / cmd_ln_int32_r(ps_search_config(ngs), "-frate");

    E_INFO("TOTAL fwdflat %.2f CPU %.3f xRT\n",
           ngs->fwdflat_perf.t_tot_cpu,
           ngs->fwdflat_perf.t_tot_cpu / n_speech);
    E_INFO("TOTAL fwdflat %.2f wall %.3f xRT\n",
           ngs->fwdflat_perf.t_tot_elapsed,
           ngs->fwdflat_perf.t_tot_elapsed / n_speech);

    /* Free single-phone words if we allocated them. */
    if (!ngs->fwdtree) {
        ngram_fwdflat_free_1ph(ngs);
    }
    ckd_free(ngs->fwdflat_wordlist);
    bitvec_free(ngs->expand_word_flag);
    ckd_free(ngs->expand_word_list);
    ckd_free(ngs->frm_wordlist);
}

int
ngram_fwdflat_reinit(ngram_search_t *ngs)
{
    /* Reallocate things that depend on the number of words. */
    int n_words;

    ckd_free(ngs->fwdflat_wordlist);
    ckd_free(ngs->expand_word_list);
    bitvec_free(ngs->expand_word_flag);
    n_words = ps_search_n_words(ngs);
    ngs->fwdflat_wordlist = ckd_calloc(n_words + 1, sizeof(*ngs->fwdflat_wordlist));
    ngs->expand_word_flag = bitvec_alloc(n_words);
    ngs->expand_word_list = ckd_calloc(n_words + 1, sizeof(*ngs->expand_word_list));
    
    /* No tree-search; take care of the expansion list and single phone words. */
    if (!ngs->fwdtree) {
        /* Free single-phone words. */
        ngram_fwdflat_free_1ph(ngs);
        /* Reallocate word_chan. */
        ckd_free(ngs->word_chan);
        ngs->word_chan = ckd_calloc(dict_size(ps_search_dict(ngs)),
                                    sizeof(*ngs->word_chan));
        /* Rebuild full expansion list from LM words. */
        ngram_fwdflat_expand_all(ngs);
        /* Allocate single phone words. */
        ngram_fwdflat_allocate_1ph(ngs);
    }
    /* Otherwise there is nothing to do since the wordlist is
     * generated anew every utterance. */
    return 0;
}

/**
 * Find all active words in backpointer table and sort by frame.
 */
static void
build_fwdflat_wordlist(ngram_search_t *ngs)
{
    int32 i, f, sf, ef, wid, nwd;
    bptbl_t *bp;
    ps_latnode_t *node, *prevnode, *nextnode;

    /* No tree-search, use statically allocated wordlist. */
    if (!ngs->fwdtree)
        return;

    memset(ngs->frm_wordlist, 0, ngs->n_frame_alloc * sizeof(*ngs->frm_wordlist));

    /* Scan the backpointer table for all active words and record
     * their exit frames. */
    for (i = 0, bp = ngs->bp_table; i < ngs->bpidx; i++, bp++) {
        sf = (bp->bp < 0) ? 0 : ngs->bp_table[bp->bp].frame + 1;
        ef = bp->frame;
        wid = bp->wid;

        /* Anything that can be transitioned to in the LM can go in
         * the word list. */
        if (!ngram_model_set_known_wid(ngs->lmset,
                                       dict_basewid(ps_search_dict(ngs), wid)))
            continue;

        /* Look for it in the wordlist. */
        for (node = ngs->frm_wordlist[sf]; node && (node->wid != wid);
             node = node->next);

        /* Update last end frame. */
        if (node)
            node->lef = ef;
        else {
            /* New node; link to head of list */
            node = listelem_malloc(ngs->latnode_alloc);
            node->wid = wid;
            node->fef = node->lef = ef;

            node->next = ngs->frm_wordlist[sf];
            ngs->frm_wordlist[sf] = node;
        }
    }

    /* Eliminate "unlikely" words, for which there are too few end points */
    for (f = 0; f < ngs->n_frame; f++) {
        prevnode = NULL;
        for (node = ngs->frm_wordlist[f]; node; node = nextnode) {
            nextnode = node->next;
            /* Word has too few endpoints */
            if ((node->lef - node->fef < ngs->min_ef_width) ||
                /* Word is </s> and doesn't actually end in last frame */
                ((node->wid == ps_search_finish_wid(ngs)) && (node->lef < ngs->n_frame - 1))) {
                if (!prevnode)
                    ngs->frm_wordlist[f] = nextnode;
                else
                    prevnode->next = nextnode;
                listelem_free(ngs->latnode_alloc, node);
            }
            else
                prevnode = node;
        }
    }

    /* Form overall wordlist for 2nd pass */
    nwd = 0;
    bitvec_clear_all(ngs->word_active, ps_search_n_words(ngs));
    for (f = 0; f < ngs->n_frame; f++) {
        for (node = ngs->frm_wordlist[f]; node; node = node->next) {
            if (!bitvec_is_set(ngs->word_active, node->wid)) {
                bitvec_set(ngs->word_active, node->wid);
                ngs->fwdflat_wordlist[nwd++] = node->wid;
            }
        }
    }
    ngs->fwdflat_wordlist[nwd] = -1;
    E_INFO("Utterance vocabulary contains %d words\n", nwd);
}

/**
 * Build HMM network for one utterance of fwdflat search.
 */
static void
build_fwdflat_chan(ngram_search_t *ngs)
{
    int32 i, wid, p;
    root_chan_t *rhmm;
    chan_t *hmm, *prevhmm;
    dict_t *dict;
    dict2pid_t *d2p;

    dict = ps_search_dict(ngs);
    d2p = ps_search_dict2pid(ngs);

    /* Build word HMMs for each word in the lattice. */
    for (i = 0; ngs->fwdflat_wordlist[i] >= 0; i++) {
        wid = ngs->fwdflat_wordlist[i];

        /* Single-phone words are permanently allocated */
        if (dict_is_single_phone(dict, wid))
            continue;

        assert(ngs->word_chan[wid] == NULL);

        /* Multiplex root HMM for first phone (one root per word, flat
         * lexicon).  diphone is irrelevant here, for the time being,
         * at least. */
        rhmm = listelem_malloc(ngs->root_chan_alloc);
        rhmm->ci2phone = dict_second_phone(dict, wid);
        rhmm->ciphone = dict_first_phone(dict, wid);
        rhmm->next = NULL;
        hmm_init(ngs->hmmctx, &rhmm->hmm, TRUE,
                 bin_mdef_pid2ssid(ps_search_acmod(ngs)->mdef, rhmm->ciphone),
                 bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef, rhmm->ciphone));

        /* HMMs for word-internal phones */
        prevhmm = NULL;
        for (p = 1; p < dict_pronlen(dict, wid) - 1; p++) {
            hmm = listelem_malloc(ngs->chan_alloc);
            hmm->ciphone = dict_pron(dict, wid, p);
            hmm->info.rc_id = (p == dict_pronlen(dict, wid) - 1) ? 0 : -1;
            hmm->next = NULL;
            hmm_init(ngs->hmmctx, &hmm->hmm, FALSE,
                     dict2pid_internal(d2p,wid,p), 
		     bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef, hmm->ciphone));

            if (prevhmm)
                prevhmm->next = hmm;
            else
                rhmm->next = hmm;

            prevhmm = hmm;
        }

        /* Right-context phones */
        ngram_search_alloc_all_rc(ngs, wid);

        /* Link in just allocated right-context phones */
        if (prevhmm)
            prevhmm->next = ngs->word_chan[wid];
        else
            rhmm->next = ngs->word_chan[wid];
        ngs->word_chan[wid] = (chan_t *) rhmm;
    }

}

void
ngram_fwdflat_start(ngram_search_t *ngs)
{
    root_chan_t *rhmm;
    int i;

    ptmr_reset(&ngs->fwdflat_perf);
    ptmr_start(&ngs->fwdflat_perf);
    build_fwdflat_wordlist(ngs);
    build_fwdflat_chan(ngs);

    ngs->bpidx = 0;
    ngs->bss_head = 0;

    for (i = 0; i < ps_search_n_words(ngs); i++)
        ngs->word_lat_idx[i] = NO_BP;

    /* Reset the permanently allocated single-phone words, since they
     * may have junk left over in them from previous searches. */
    for (i = 0; i < ngs->n_1ph_words; i++) {
        int32 w = ngs->single_phone_wid[i];
        rhmm = (root_chan_t *) ngs->word_chan[w];
        hmm_clear(&rhmm->hmm);
    }

    /* Start search with <s>; word_chan[<s>] is permanently allocated */
    rhmm = (root_chan_t *) ngs->word_chan[ps_search_start_wid(ngs)];
    hmm_enter(&rhmm->hmm, 0, NO_BP, 0);
    ngs->active_word_list[0][0] = ps_search_start_wid(ngs);
    ngs->n_active_word[0] = 1;

    ngs->best_score = 0;
    ngs->renormalized = FALSE;

    for (i = 0; i < ps_search_n_words(ngs); i++)
        ngs->last_ltrans[i].sf = -1;

    if (!ngs->fwdtree)
        ngs->n_frame = 0;

    ngs->st.n_fwdflat_chan = 0;
    ngs->st.n_fwdflat_words = 0;
    ngs->st.n_fwdflat_word_transition = 0;
    ngs->st.n_senone_active_utt = 0;
}

static void
compute_fwdflat_sen_active(ngram_search_t *ngs, int frame_idx)
{
    int32 i, w;
    int32 *awl;
    root_chan_t *rhmm;
    chan_t *hmm;

    acmod_clear_active(ps_search_acmod(ngs));

    i = ngs->n_active_word[frame_idx & 0x1];
    awl = ngs->active_word_list[frame_idx & 0x1];

    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (root_chan_t *)ngs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            acmod_activate_hmm(ps_search_acmod(ngs), &rhmm->hmm);
        }

        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) == frame_idx) {
                acmod_activate_hmm(ps_search_acmod(ngs), &hmm->hmm);
            }
        }
    }
}

static void
fwdflat_eval_chan(ngram_search_t *ngs, int frame_idx)
{
    int32 i, w, bestscore;
    int32 *awl;
    root_chan_t *rhmm;
    chan_t *hmm;

    i = ngs->n_active_word[frame_idx & 0x1];
    awl = ngs->active_word_list[frame_idx & 0x1];
    bestscore = WORST_SCORE;

    ngs->st.n_fwdflat_words += i;

    /* Scan all active words. */
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (root_chan_t *) ngs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            int32 score = chan_v_eval(rhmm);
            if ((score BETTER_THAN bestscore) && (w != ps_search_finish_wid(ngs)))
                bestscore = score;
            ngs->st.n_fwdflat_chan++;
        }

        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) == frame_idx) {
                int32 score = chan_v_eval(hmm);
                if (score BETTER_THAN bestscore)
                    bestscore = score;
                ngs->st.n_fwdflat_chan++;
            }
        }
    }

    ngs->best_score = bestscore;
}

static void
fwdflat_prune_chan(ngram_search_t *ngs, int frame_idx)
{
    int32 i, cf, nf, w, pip, newscore, thresh, wordthresh;
    int32 *awl;
    root_chan_t *rhmm;
    chan_t *hmm, *nexthmm;

    cf = frame_idx;
    nf = cf + 1;
    i = ngs->n_active_word[cf & 0x1];
    awl = ngs->active_word_list[cf & 0x1];
    bitvec_clear_all(ngs->word_active, ps_search_n_words(ngs));

    thresh = ngs->best_score + ngs->fwdflatbeam;
    wordthresh = ngs->best_score + ngs->fwdflatwbeam;
    pip = ngs->pip;
    E_DEBUG(3,("frame %d thresh %d wordthresh %d\n", frame_idx, thresh, wordthresh));

    /* Scan all active words. */
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (root_chan_t *) ngs->word_chan[w];
        /* Propagate active root channels */
        if (hmm_frame(&rhmm->hmm) == cf
            && hmm_bestscore(&rhmm->hmm) BETTER_THAN thresh) {
            hmm_frame(&rhmm->hmm) = nf;
            bitvec_set(ngs->word_active, w);

            /* Transitions out of root channel */
            newscore = hmm_out_score(&rhmm->hmm);
            if (rhmm->next) {
                assert(!dict_is_single_phone(ps_search_dict(ngs), w));

                newscore += pip;
                if (newscore BETTER_THAN thresh) {
                    hmm = rhmm->next;
                    /* Enter all right context phones */
                    if (hmm->info.rc_id >= 0) {
                        for (; hmm; hmm = hmm->next) {
                            if ((hmm_frame(&hmm->hmm) < cf)
                                || (newscore BETTER_THAN hmm_in_score(&hmm->hmm))) {
                                hmm_enter(&hmm->hmm, newscore,
                                          hmm_out_history(&rhmm->hmm), nf);
                            }
                        }
                    }
                    /* Just a normal word internal phone */
                    else {
                        if ((hmm_frame(&hmm->hmm) < cf)
                            || (newscore BETTER_THAN hmm_in_score(&hmm->hmm))) {
                                hmm_enter(&hmm->hmm, newscore,
                                          hmm_out_history(&rhmm->hmm), nf);
                        }
                    }
                }
            }
            else {
                assert(dict_is_single_phone(ps_search_dict(ngs), w));

                /* Word exit for single-phone words (where did their
                 * whmms come from?) (either from
                 * ngram_search_fwdtree, or from
                 * ngram_fwdflat_allocate_1ph(), that's where) */
                if (newscore BETTER_THAN wordthresh) {
                    ngram_search_save_bp(ngs, cf, w, newscore,
                                         hmm_out_history(&rhmm->hmm), 0);
                }
            }
        }

        /* Transitions out of non-root channels. */
        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) >= cf) {
                /* Propagate forward HMMs inside the beam. */
                if (hmm_bestscore(&hmm->hmm) BETTER_THAN thresh) {
                    hmm_frame(&hmm->hmm) = nf;
                    bitvec_set(ngs->word_active, w);

                    newscore = hmm_out_score(&hmm->hmm);
                    /* Word-internal phones */
                    if (hmm->info.rc_id < 0) {
                        newscore += pip;
                        if (newscore BETTER_THAN thresh) {
                            nexthmm = hmm->next;
                            /* Enter all right-context phones. */
                            if (nexthmm->info.rc_id >= 0) {
                                 for (; nexthmm; nexthmm = nexthmm->next) {
                                    if ((hmm_frame(&nexthmm->hmm) < cf)
                                        || (newscore BETTER_THAN
                                            hmm_in_score(&nexthmm->hmm))) {
                                        hmm_enter(&nexthmm->hmm,
                                                  newscore,
                                                  hmm_out_history(&hmm->hmm),
                                                  nf);
                                    }
                                }
                            }
                            /* Enter single word-internal phone. */
                            else {
                                if ((hmm_frame(&nexthmm->hmm) < cf)
                                    || (newscore BETTER_THAN
                                        hmm_in_score(&nexthmm->hmm))) {
                                    hmm_enter(&nexthmm->hmm, newscore,
                                              hmm_out_history(&hmm->hmm), nf);
                                }
                            }
                        }
                    }
                    /* Right-context phones - apply word beam and exit. */
                    else {
                        if (newscore BETTER_THAN wordthresh) {
                            ngram_search_save_bp(ngs, cf, w, newscore,
                                                 hmm_out_history(&hmm->hmm),
                                                 hmm->info.rc_id);
                        }
                    }
                }
                /* Zero out inactive HMMs. */
                else if (hmm_frame(&hmm->hmm) != nf) {
                    hmm_clear_scores(&hmm->hmm);
                }
            }
        }
    }
}

static void
get_expand_wordlist(ngram_search_t *ngs, int32 frm, int32 win)
{
    int32 f, sf, ef;
    ps_latnode_t *node;

    if (!ngs->fwdtree) {
        ngs->st.n_fwdflat_word_transition += ngs->n_expand_words;
        return;
    }

    sf = frm - win;
    if (sf < 0)
        sf = 0;
    ef = frm + win;
    if (ef > ngs->n_frame)
        ef = ngs->n_frame;

    bitvec_clear_all(ngs->expand_word_flag, ps_search_n_words(ngs));
    ngs->n_expand_words = 0;

    for (f = sf; f < ef; f++) {
        for (node = ngs->frm_wordlist[f]; node; node = node->next) {
            if (!bitvec_is_set(ngs->expand_word_flag, node->wid)) {
                ngs->expand_word_list[ngs->n_expand_words++] = node->wid;
                bitvec_set(ngs->expand_word_flag, node->wid);
            }
        }
    }
    ngs->expand_word_list[ngs->n_expand_words] = -1;
    ngs->st.n_fwdflat_word_transition += ngs->n_expand_words;
}

static void
fwdflat_word_transition(ngram_search_t *ngs, int frame_idx)
{
    int32 cf, nf, b, thresh, pip, i, w, newscore;
    int32 best_silrc_score = 0, best_silrc_bp = 0;      /* FIXME: good defaults? */
    bptbl_t *bp;
    int32 *rcss;
    root_chan_t *rhmm;
    int32 *awl;
    float32 lwf;
    dict_t *dict = ps_search_dict(ngs);
    dict2pid_t *d2p = ps_search_dict2pid(ngs);

    cf = frame_idx;
    nf = cf + 1;
    thresh = ngs->best_score + ngs->fwdflatbeam;
    pip = ngs->pip;
    best_silrc_score = WORST_SCORE;
    lwf = ngs->fwdflat_fwdtree_lw_ratio;

    /* Search for all words starting within a window of this frame.
     * These are the successors for words exiting now. */
    get_expand_wordlist(ngs, cf, ngs->max_sf_win);

    /* Scan words exited in current frame */
    for (b = ngs->bp_table_idx[cf]; b < ngs->bpidx; b++) {
        xwdssid_t *rssid;
        int32 silscore;

        bp = ngs->bp_table + b;
        ngs->word_lat_idx[bp->wid] = NO_BP;

        if (bp->wid == ps_search_finish_wid(ngs))
            continue;

        /* DICT2PID location */
        /* Get the mapping from right context phone ID to index in the
         * right context table and the bscore_stack. */
        rcss = ngs->bscore_stack + bp->s_idx;
        if (bp->last2_phone == -1)
            rssid = NULL;
        else
            rssid = dict2pid_rssid(d2p, bp->last_phone, bp->last2_phone);

        /* Transition to all successor words. */
        for (i = 0; ngs->expand_word_list[i] >= 0; i++) {
            int32 n_used;

            w = ngs->expand_word_list[i];

            /* Get the exit score we recorded in save_bwd_ptr(), or
             * something approximating it. */
            if (rssid)
                newscore = rcss[rssid->cimap[dict_first_phone(dict, w)]];
            else
                newscore = bp->score;
            if (newscore == WORST_SCORE)
                continue;
            /* FIXME: Floating point... */
            newscore += lwf
                * (ngram_tg_score(ngs->lmset,
                                  dict_basewid(dict, w),
                                  bp->real_wid,
                                  bp->prev_real_wid,
                                  &n_used) >> SENSCR_SHIFT);
            newscore += pip;

            /* Enter the next word */
            if (newscore BETTER_THAN thresh) {
                rhmm = (root_chan_t *) ngs->word_chan[w];
                if ((hmm_frame(&rhmm->hmm) < cf)
                    || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                    hmm_enter(&rhmm->hmm, newscore, b, nf);
                    /* DICT2PID: This is where mpx ssids get introduced. */
                    /* Look up the ssid to use when entering this mpx triphone. */
                    hmm_mpx_ssid(&rhmm->hmm, 0) =
                        dict2pid_ldiph_lc(d2p, rhmm->ciphone, rhmm->ci2phone,
                                          dict_last_phone(dict, bp->wid));
                    assert(IS_S3SSID(hmm_mpx_ssid(&rhmm->hmm, 0)));
                    E_DEBUG(6,("ssid %d(%d,%d) = %d\n",
                               rhmm->ciphone, dict_last_phone(dict, bp->wid), rhmm->ci2phone,
                               hmm_mpx_ssid(&rhmm->hmm, 0)));
                    bitvec_set(ngs->word_active, w);
                }
            }
        }

        /* Get the best exit into silence. */
        if (rssid)
            silscore = rcss[rssid->cimap[ps_search_acmod(ngs)->mdef->sil]];
        else
            silscore = bp->score;
        if (silscore BETTER_THAN best_silrc_score) {
            best_silrc_score = silscore;
            best_silrc_bp = b;
        }
    }

    /* Transition to <sil> */
    newscore = best_silrc_score + ngs->silpen + pip;
    if ((newscore BETTER_THAN thresh) && (newscore BETTER_THAN WORST_SCORE)) {
        w = ps_search_silence_wid(ngs);
        rhmm = (root_chan_t *) ngs->word_chan[w];
        if ((hmm_frame(&rhmm->hmm) < cf)
            || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
            hmm_enter(&rhmm->hmm, newscore,
                      best_silrc_bp, nf);
            bitvec_set(ngs->word_active, w);
        }
    }
    /* Transition to noise words */
    newscore = best_silrc_score + ngs->fillpen + pip;
    if ((newscore BETTER_THAN thresh) && (newscore BETTER_THAN WORST_SCORE)) {
        for (w = ps_search_silence_wid(ngs) + 1; w < ps_search_n_words(ngs); w++) {
            rhmm = (root_chan_t *) ngs->word_chan[w];
            /* Noise words that aren't a single phone will have NULL here. */
            if (rhmm == NULL)
                continue;
            if ((hmm_frame(&rhmm->hmm) < cf)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm, newscore,
                          best_silrc_bp, nf);
                bitvec_set(ngs->word_active, w);
            }
        }
    }

    /* Reset initial channels of words that have become inactive even after word trans. */
    i = ngs->n_active_word[cf & 0x1];
    awl = ngs->active_word_list[cf & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (root_chan_t *) ngs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == cf) {
            hmm_clear_scores(&rhmm->hmm);
        }
    }
}

static void
fwdflat_renormalize_scores(ngram_search_t *ngs, int frame_idx, int32 norm)
{
    root_chan_t *rhmm;
    chan_t *hmm;
    int32 i, cf, w, *awl;

    cf = frame_idx;

    /* Renormalize individual word channels */
    i = ngs->n_active_word[cf & 0x1];
    awl = ngs->active_word_list[cf & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (root_chan_t *) ngs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == cf) {
            hmm_normalize(&rhmm->hmm, norm);
        }
        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) == cf) {
                hmm_normalize(&hmm->hmm, norm);
            }
        }
    }

    ngs->renormalized = TRUE;
}

int
ngram_fwdflat_search(ngram_search_t *ngs, int frame_idx)
{
    int16 const *senscr;
    int32 nf, i, j;
    int32 *nawl;

    /* Activate our HMMs for the current frame if need be. */
    if (!ps_search_acmod(ngs)->compallsen)
        compute_fwdflat_sen_active(ngs, frame_idx);

    /* Compute GMM scores for the current frame. */
    senscr = acmod_score(ps_search_acmod(ngs), &frame_idx);
    ngs->st.n_senone_active_utt += ps_search_acmod(ngs)->n_senone_active;

    /* Mark backpointer table for current frame. */
    ngram_search_mark_bptable(ngs, frame_idx);

    /* If the best score is equal to or worse than WORST_SCORE,
     * recognition has failed, don't bother to keep trying. */
    if (ngs->best_score == WORST_SCORE || ngs->best_score WORSE_THAN WORST_SCORE)
        return 0;
    /* Renormalize if necessary */
    if (ngs->best_score + (2 * ngs->beam) WORSE_THAN WORST_SCORE) {
        E_INFO("Renormalizing Scores at frame %d, best score %d\n",
               frame_idx, ngs->best_score);
        fwdflat_renormalize_scores(ngs, frame_idx, ngs->best_score);
    }

    ngs->best_score = WORST_SCORE;
    hmm_context_set_senscore(ngs->hmmctx, senscr);

    /* Evaluate HMMs */
    fwdflat_eval_chan(ngs, frame_idx);
    /* Prune HMMs and do phone transitions. */
    fwdflat_prune_chan(ngs, frame_idx);
    /* Do word transitions. */
    fwdflat_word_transition(ngs, frame_idx);

    /* Create next active word list */
    nf = frame_idx + 1;
    nawl = ngs->active_word_list[nf & 0x1];
    for (i = 0, j = 0; ngs->fwdflat_wordlist[i] >= 0; i++) {
        if (bitvec_is_set(ngs->word_active, ngs->fwdflat_wordlist[i])) {
            *(nawl++) = ngs->fwdflat_wordlist[i];
            j++;
        }
    }
    for (i = ps_search_start_wid(ngs); i < ps_search_n_words(ngs); i++) {
        if (bitvec_is_set(ngs->word_active, i)) {
            *(nawl++) = i;
            j++;
        }
    }
    if (!ngs->fwdtree)
        ++ngs->n_frame;
    ngs->n_active_word[nf & 0x1] = j;

    /* Return the number of frames processed. */
    return 1;
}

/**
 * Destroy wordlist from the current utterance.
 */
static void
destroy_fwdflat_wordlist(ngram_search_t *ngs)
{
    ps_latnode_t *node, *tnode;
    int32 f;

    if (!ngs->fwdtree)
        return;

    for (f = 0; f < ngs->n_frame; f++) {
        for (node = ngs->frm_wordlist[f]; node; node = tnode) {
            tnode = node->next;
            listelem_free(ngs->latnode_alloc, node);
        }
    }
}

/**
 * Free HMM network for one utterance of fwdflat search.
 */
static void
destroy_fwdflat_chan(ngram_search_t *ngs)
{
    int32 i, wid;

    for (i = 0; ngs->fwdflat_wordlist[i] >= 0; i++) {
        root_chan_t *rhmm;
        chan_t *thmm;
        wid = ngs->fwdflat_wordlist[i];
        if (dict_is_single_phone(ps_search_dict(ngs),wid))
            continue;
        assert(ngs->word_chan[wid] != NULL);

        /* The first HMM in ngs->word_chan[wid] was allocated with
         * ngs->root_chan_alloc, but this will attempt to free it
         * using ngs->chan_alloc, which will not work.  Therefore we
         * free it manually and move the list forward before handing
         * it off. */
        rhmm = (root_chan_t *)ngs->word_chan[wid];
        thmm = rhmm->next;
        listelem_free(ngs->root_chan_alloc, rhmm);
        ngs->word_chan[wid] = thmm;
        ngram_search_free_all_rc(ngs, wid);
    }
}

void
ngram_fwdflat_finish(ngram_search_t *ngs)
{
    int32 cf;

    destroy_fwdflat_chan(ngs);
    destroy_fwdflat_wordlist(ngs);
    bitvec_clear_all(ngs->word_active, ps_search_n_words(ngs));

    /* This is the number of frames processed. */
    cf = ps_search_acmod(ngs)->output_frame;
    /* Add a mark in the backpointer table for one past the final frame. */
    ngram_search_mark_bptable(ngs, cf);

    ptmr_stop(&ngs->fwdflat_perf);
    /* Print out some statistics. */
    if (cf > 0) {
        double n_speech = (double)(cf + 1)
            / cmd_ln_int32_r(ps_search_config(ngs), "-frate");
        E_INFO("%8d words recognized (%d/fr)\n",
               ngs->bpidx, (ngs->bpidx + (cf >> 1)) / (cf + 1));
        E_INFO("%8d senones evaluated (%d/fr)\n", ngs->st.n_senone_active_utt,
               (ngs->st.n_senone_active_utt + (cf >> 1)) / (cf + 1));
        E_INFO("%8d channels searched (%d/fr)\n",
               ngs->st.n_fwdflat_chan, ngs->st.n_fwdflat_chan / (cf + 1));
        E_INFO("%8d words searched (%d/fr)\n",
               ngs->st.n_fwdflat_words, ngs->st.n_fwdflat_words / (cf + 1));
        E_INFO("%8d word transitions (%d/fr)\n",
               ngs->st.n_fwdflat_word_transition,
               ngs->st.n_fwdflat_word_transition / (cf + 1));
        E_INFO("fwdflat %.2f CPU %.3f xRT\n",
               ngs->fwdflat_perf.t_cpu,
               ngs->fwdflat_perf.t_cpu / n_speech);
        E_INFO("fwdflat %.2f wall %.3f xRT\n",
               ngs->fwdflat_perf.t_elapsed,
               ngs->fwdflat_perf.t_elapsed / n_speech);
    }
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ngram_search_fwdflat.h Flat lexicon based Viterbi search.
 */

#ifndef __NGRAM_SEARCH_FWDFLAT_H__
#define __NGRAM_SEARCH_FWDFLAT_H__

/* SphinxBase headers. */

/* Local headers. */
#include "ngram_search.h"

/**
 * Initialize N-Gram search for fwdflat decoding.
 */
void ngram_fwdflat_init(ngram_search_t *ngs);

/**
 * Release memory associated with fwdflat decoding.
 */
void ngram_fwdflat_deinit(ngram_search_t *ngs);

/**
 * Rebuild search structures for updated language models.
 */
int ngram_fwdflat_reinit(ngram_search_t *ngs);

/**
 * Start fwdflat decoding for an utterance.
 */
void ngram_fwdflat_start(ngram_search_t *ngs);

/**
 * Search one frame forward in an utterance.
 */
int ngram_fwdflat_search(ngram_search_t *ngs, int frame_idx);

/**
 * Finish fwdflat decoding for an utterance.
 */
void ngram_fwdflat_finish(ngram_search_t *ngs);


#endif /* __NGRAM_SEARCH_FWDFLAT_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ngram_search_fwdtree.c Lexicon tree search.
 */

/* System headers. */
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "ngram_search_fwdtree.h"
#include "phone_loop_search.h"

/* Turn this on to dump channels for debugging */
#define __CHAN_DUMP__		0
#if __CHAN_DUMP__
#define chan_v_eval(chan) hmm_dump_vit_eval(&(chan)->hmm, stderr)
#else
#define chan_v_eval(chan) hmm_vit_eval(&(chan)->hmm)
#endif

/*
 * Allocate that part of the search channel tree structure that is independent of the
 * LM in use.
 */
static void
init_search_tree(ngram_search_t *ngs)
{
    int32 w, ndiph, i, n_words, n_ci;
    dict_t *dict = ps_search_dict(ngs);
    bitvec_t *dimap;

    n_words = ps_search_n_words(ngs);
    ngs->homophone_set = ckd_calloc(n_words, sizeof(*ngs->homophone_set));

    /* Find #single phone words, and #unique first diphones (#root channels) in dict. */
    ndiph = 0;
    ngs->n_1ph_words = 0;
    n_ci = bin_mdef_n_ciphone(ps_search_acmod(ngs)->mdef);
    /* Allocate a bitvector with flags for each possible diphone. */
    dimap = bitvec_alloc(n_ci * n_ci);
    for (w = 0; w < n_words; w++) {
        if (!dict_real_word(dict, w))
            continue;
        if (dict_is_single_phone(dict, w))
            ++ngs->n_1ph_words;
        else {
            int ph0, ph1;
            ph0 = dict_first_phone(dict, w);
            ph1 = dict_second_phone(dict, w);
            /* Increment ndiph the first time we see a diphone. */
            if (bitvec_is_clear(dimap, ph0 * n_ci + ph1)) {
                bitvec_set(dimap, ph0 * n_ci + ph1);
                ++ndiph;
            }
        }
    }
    E_INFO("%d unique initial diphones\n", ndiph);
    bitvec_free(dimap);

    /* Add remaining dict words (</s>, <s>, <sil>, noise words) to single-phone words */
    ngs->n_1ph_words += dict_num_fillers(dict) + 2;
    ngs->n_root_chan_alloc = ndiph + 1;
    /* Verify that these are all *actually* single-phone words,
     * otherwise really bad things will happen to us. */
    for (w = 0; w < n_words; ++w) {
        if (dict_real_word(dict, w))
            continue;
        if (!dict_is_single_phone(dict, w)) {
            E_WARN("Filler word %d = %s has more than one phone, ignoring it.\n",
                   w, dict_wordstr(dict, w));
            --ngs->n_1ph_words;
        }
    }

    /* Allocate and initialize root channels */
    ngs->root_chan =
        ckd_calloc(ngs->n_root_chan_alloc, sizeof(*ngs->root_chan));
    for (i = 0; i < ngs->n_root_chan_alloc; i++) {
        hmm_init(ngs->hmmctx, &ngs->root_chan[i].hmm, TRUE, -1, -1);
        ngs->root_chan[i].penult_phn_wid = -1;
        ngs->root_chan[i].next = NULL;
    }

    /* Permanently allocate and initialize channels for single-phone
     * words (1/word). */
    ngs->rhmm_1ph = ckd_calloc(ngs->n_1ph_words, sizeof(*ngs->rhmm_1ph));
    i = 0;
    for (w = 0; w < n_words; w++) {
        if (!dict_is_single_phone(dict, w))
            continue;
        /* Use SIL as right context for these. */
        ngs->rhmm_1ph[i].ci2phone = bin_mdef_silphone(ps_search_acmod(ngs)->mdef);
        ngs->rhmm_1ph[i].ciphone = dict_first_phone(dict, w);
        hmm_init(ngs->hmmctx, &ngs->rhmm_1ph[i].hmm, TRUE,
                 bin_mdef_pid2ssid(ps_search_acmod(ngs)->mdef, ngs->rhmm_1ph[i].ciphone),
                 bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef, ngs->rhmm_1ph[i].ciphone));
        ngs->rhmm_1ph[i].next = NULL;

        ngs->word_chan[w] = (chan_t *) &(ngs->rhmm_1ph[i]);
        i++;
    }

    ngs->single_phone_wid = ckd_calloc(ngs->n_1ph_words,
                                       sizeof(*ngs->single_phone_wid));
    E_INFO("%d root, %d non-root channels, %d single-phone words\n",
           ngs->n_root_chan, ngs->n_nonroot_chan, ngs->n_1ph_words);
}

/*
 * One-time initialization of internal channels in HMM tree.
 */
static void
init_nonroot_chan(ngram_search_t *ngs, chan_t * hmm, int32 ph, int32 ci, int32 tmatid)
{
    hmm->next = NULL;
    hmm->alt = NULL;
    hmm->info.penult_phn_wid = -1;
    hmm->ciphone = ci;
    hmm_init(ngs->hmmctx, &hmm->hmm, FALSE, ph, tmatid);
}

/*
 * Allocate and initialize search channel-tree structure.
 * At this point, all the root-channels have been allocated and partly initialized
 * (as per init_search_tree()), and channels for all the single-phone words have been
 * allocated and initialized.  None of the interior channels of search-trees have
 * been allocated.
 * This routine may be called on every utterance, after reinit_search_tree() clears
 * the search tree created for the previous utterance.  Meant for reconfiguring the
 * search tree to suit the currently active LM.
 */
static void
create_search_tree(ngram_search_t *ngs)
{
    chan_t *hmm;
    root_chan_t *rhmm;
    int32 w, i, j, p, ph, tmatid;
    int32 n_words;
    dict_t *dict = ps_search_dict(ngs);
    dict2pid_t *d2p = ps_search_dict2pid(ngs);

    n_words = ps_search_n_words(ngs);

    E_INFO("Creating search tree\n");

    for (w = 0; w < n_words; w++)
        ngs->homophone_set[w] = -1;

    E_INFO("before: %d root, %d non-root channels, %d single-phone words\n",
           ngs->n_root_chan, ngs->n_nonroot_chan, ngs->n_1ph_words);

    ngs->n_1ph_LMwords = 0;
    ngs->n_root_chan = 0;
    ngs->n_nonroot_chan = 0;

    for (w = 0; w < n_words; w++) {
        int ciphone, ci2phone;

        /* Ignore dictionary words not in LM */
        if (!ngram_model_set_known_wid(ngs->lmset, dict_basewid(dict, w)))
            continue;

        /* Handle single-phone words individually; not in channel tree */
        if (dict_is_single_phone(dict, w)) {
            E_DEBUG(1,("single_phone_wid[%d] = %s\n",
                       ngs->n_1ph_LMwords, dict_wordstr(dict, w)));
            ngs->single_phone_wid[ngs->n_1ph_LMwords++] = w;
            continue;
        }

        /* Find a root channel matching the initial diphone, or
         * allocate one if not found. */
        ciphone = dict_first_phone(dict, w);
        ci2phone = dict_second_phone(dict, w);
        for (i = 0; i < ngs->n_root_chan; ++i) {
            if (ngs->root_chan[i].ciphone == ciphone
                && ngs->root_chan[i].ci2phone == ci2phone)
                break;
        }
        if (i == ngs->n_root_chan) {
            rhmm = &(ngs->root_chan[ngs->n_root_chan]);
            rhmm->hmm.tmatid = bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef, ciphone);
            /* Begin with CI phone?  Not sure this makes a difference... */
            hmm_mpx_ssid(&rhmm->hmm, 0) =
                bin_mdef_pid2ssid(ps_search_acmod(ngs)->mdef, ciphone);
            rhmm->ciphone = ciphone;
            rhmm->ci2phone = ci2phone;
            ngs->n_root_chan++;
        }
        else
            rhmm = &(ngs->root_chan[i]);

        E_DEBUG(3,("word %s rhmm %d\n", dict_wordstr(dict, w), rhmm - ngs->root_chan));
        /* Now, rhmm = root channel for w.  Go on to remaining phones */
        if (dict_pronlen(dict, w) == 2) {
            /* Next phone is the last; not kept in tree; add w to penult_phn_wid set */
            if ((j = rhmm->penult_phn_wid) < 0)
                rhmm->penult_phn_wid = w;
            else {
                for (; ngs->homophone_set[j] >= 0; j = ngs->homophone_set[j]);
                ngs->homophone_set[j] = w;
            }
        }
        else {
            /* Add remaining phones, except the last, to tree */
            ph = dict2pid_internal(d2p, w, 1);
            tmatid = bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef, dict_pron(dict, w, 1));
            hmm = rhmm->next;
            if (hmm == NULL) {
                rhmm->next = hmm = listelem_malloc(ngs->chan_alloc);
                init_nonroot_chan(ngs, hmm, ph, dict_pron(dict, w, 1), tmatid);
                ngs->n_nonroot_chan++;
            }
            else {
                chan_t *prev_hmm = NULL;

                for (; hmm && (hmm_nonmpx_ssid(&hmm->hmm) != ph); hmm = hmm->alt)
                    prev_hmm = hmm;
                if (!hmm) {     /* thanks, rkm! */
                    prev_hmm->alt = hmm = listelem_malloc(ngs->chan_alloc);
                    init_nonroot_chan(ngs, hmm, ph, dict_pron(dict, w, 1), tmatid);
                    ngs->n_nonroot_chan++;
                }
            }
            E_DEBUG(3,("phone %s = %d\n",
                       bin_mdef_ciphone_str(ps_search_acmod(ngs)->mdef,
                                            dict_second_phone(dict, w)), ph));
            for (p = 2; p < dict_pronlen(dict, w) - 1; p++) {
                ph = dict2pid_internal(d2p, w, p);
                tmatid = bin_mdef_pid2tmatid(ps_search_acmod(ngs)->mdef, dict_pron(dict, w, p));
                if (!hmm->next) {
                    hmm->next = listelem_malloc(ngs->chan_alloc);
                    hmm = hmm->next;
                    init_nonroot_chan(ngs, hmm, ph, dict_pron(dict, w, p), tmatid);
                    ngs->n_nonroot_chan++;
                }
                else {
                    chan_t *prev_hmm = NULL;

                    for (hmm = hmm->next; hmm && (hmm_nonmpx_ssid(&hmm->hmm) != ph);
                         hmm = hmm->alt)
                        prev_hmm = hmm;
                    if (!hmm) { /* thanks, rkm! */
                        prev_hmm->alt = hmm = listelem_malloc(ngs->chan_alloc);
                        init_nonroot_chan(ngs, hmm, ph, dict_pron(dict, w, p), tmatid);
                        ngs->n_nonroot_chan++;
                    }
                }
                E_DEBUG(3,("phone %s = %d\n",
                           bin_mdef_ciphone_str(ps_search_acmod(ngs)->mdef,
                                                dict_pron(dict, w, p)), ph));
            }

            /* All but last phone of w in tree; add w to hmm->info.penult_phn_wid set */
            if ((j = hmm->info.penult_phn_wid) < 0)
                hmm->info.penult_phn_wid = w;
            else {
                for (; ngs->homophone_set[j] >= 0; j = ngs->homophone_set[j]);
                ngs->homophone_set[j] = w;
            }
        }
    }

    ngs->n_1ph_words = ngs->n_1ph_LMwords;

    /* Add filler words to the array of 1ph words. */
    for (w = 0; w < n_words; ++w) {
        /* Skip anything that doesn't actually have a single phone. */
        if (!dict_is_single_phone(dict, w))
            continue;
        /* Also skip "real words" and things that are in the LM. */
        if (dict_real_word(dict, w))
            continue;
        if (ngram_model_set_known_wid(ngs->lmset, dict_basewid(dict, w)))
            continue;
        E_DEBUG(1,("single_phone_wid[%d] = %s\n",
                   ngs->n_1ph_words, dict_wordstr(dict, w)));
        ngs->single_phone_wid[ngs->n_1ph_words++] = w;
    }

    if (ngs->n_nonroot_chan >= ngs->max_nonroot_chan) {
        /* Give some room for channels for new words added dynamically at run time */
        ngs->max_nonroot_chan = ngs->n_nonroot_chan + 128;
        E_INFO("after: max nonroot chan increased to %d\n", ngs->max_nonroot_chan);

        /* Free old active channel list array if any and allocate new one */
        if (ngs->active_chan_list)
            ckd_free_2d(ngs->active_chan_list);
        ngs->active_chan_list = ckd_calloc_2d(2, ngs->max_nonroot_chan,
                                              sizeof(**ngs->active_chan_list));
    }

    if (!ngs->n_root_chan)
	E_ERROR("No word from the language model has pronunciation in the dictionary\n");

    E_INFO("after: %d root, %d non-root channels, %d single-phone words\n",
           ngs->n_root_chan, ngs->n_nonroot_chan, ngs->n_1ph_words);
}

static void
reinit_search_subtree(ngram_search_t *ngs, chan_t * hmm)
{
    chan_t *child, *sibling;

    /* First free all children under hmm */
    for (child = hmm->next; child; child = sibling) {
        sibling = child->alt;
        reinit_search_subtree(ngs, child);
    }

    /* Now free hmm */
    hmm_deinit(&hmm->hmm);
    listelem_free(ngs->chan_alloc, hmm);
}

/*
 * Delete search tree by freeing all interior channels within search tree and
 * restoring root channel state to the init state (i.e., just after init_search_tree()).
 */
static void
reinit_search_tree(ngram_search_t *ngs)
{
    int32 i;
    chan_t *hmm, *sibling;

    for (i = 0; i < ngs->n_root_chan; i++) {
        hmm = ngs->root_chan[i].next;

        while (hmm) {
            sibling = hmm->alt;
            reinit_search_subtree(ngs, hmm);
            hmm = sibling;
        }

        ngs->root_chan[i].penult_phn_wid = -1;
        ngs->root_chan[i].next = NULL;
    }
    ngs->n_nonroot_chan = 0;
}

void
ngram_fwdtree_init(ngram_search_t *ngs)
{
    /* Allocate bestbp_rc, lastphn_cand, last_ltrans */
    ngs->bestbp_rc = ckd_calloc(bin_mdef_n_ciphone(ps_search_acmod(ngs)->mdef),
                                sizeof(*ngs->bestbp_rc));
    ngs->lastphn_cand = ckd_calloc(ps_search_n_words(ngs),
                                   sizeof(*ngs->lastphn_cand));
    init_search_tree(ngs);
    create_search_tree(ngs);
}

static void
deinit_search_tree(ngram_search_t *ngs)
{
    int i, w, n_words;

    n_words = ps_search_n_words(ngs);
    for (i = 0; i < ngs->n_root_chan_alloc; i++) {
        hmm_deinit(&ngs->root_chan[i].hmm);
    }
    if (ngs->rhmm_1ph) {
        for (i = w = 0; w < n_words; ++w) {
            if (!dict_is_single_phone(ps_search_dict(ngs), w))
                continue;
            hmm_deinit(&ngs->rhmm_1ph[i].hmm);
            ++i;
        }
        ckd_free(ngs->rhmm_1ph);
        ngs->rhmm_1ph = NULL;
    }
    ngs->n_root_chan = 0;
    ngs->n_root_chan_alloc = 0;
    ckd_free(ngs->root_chan);
    ngs->root_chan = NULL;
    ckd_free(ngs->single_phone_wid);
    ngs->single_phone_wid = NULL;
    ckd_free(ngs->homophone_set);
    ngs->homophone_set = NULL;
}

void
ngram_fwdtree_deinit(ngram_search_t *ngs)
{
    double n_speech = (double)ngs->n_tot_frame
            / cmd_ln_int32_r(ps_search_config(ngs), "-frate");

    E_INFO("TOTAL fwdtree %.2f CPU %.3f xRT\n",
           ngs->fwdtree_perf.t_tot_cpu,
           ngs->fwdtree_perf.t_tot_cpu / n_speech);
    E_INFO("TOTAL fwdtree %.2f wall %.3f xRT\n",
           ngs->fwdtree_perf.t_tot_elapsed,
           ngs->fwdtree_perf.t_tot_elapsed / n_speech);

    /* Reset non-root channels. */
    reinit_search_tree(ngs);
    /* Free the search tree. */
    deinit_search_tree(ngs);
    /* Free other stuff. */
    ngs->max_nonroot_chan = 0;
    ckd_free_2d(ngs->active_chan_list);
    ngs->active_chan_list = NULL;
    ckd_free(ngs->cand_sf);
    ngs->cand_sf = NULL;
    ckd_free(ngs->bestbp_rc);
    ngs->bestbp_rc = NULL;
    ckd_free(ngs->lastphn_cand);
    ngs->lastphn_cand = NULL;
}

int
ngram_fwdtree_reinit(ngram_search_t *ngs)
{
    /* Reset non-root channels. */
    reinit_search_tree(ngs);
    /* Free the search tree. */
    deinit_search_tree(ngs);
    /* Reallocate things that depend on the number of words. */
    ckd_free(ngs->lastphn_cand);
    ngs->lastphn_cand = ckd_calloc(ps_search_n_words(ngs),
                                   sizeof(*ngs->lastphn_cand));
    ckd_free(ngs->word_chan);
    ngs->word_chan = ckd_calloc(ps_search_n_words(ngs),
                                sizeof(*ngs->word_chan));
    /* Rebuild the search tree. */
    init_search_tree(ngs);
    create_search_tree(ngs);
    return 0;
}

void
ngram_fwdtree_start(ngram_search_t *ngs)
{
    ps_search_t *base = (ps_search_t *)ngs;
    int32 i, w, n_words;
    root_chan_t *rhmm;

    n_words = ps_search_n_words(ngs);

    /* Reset utterance statistics. */
    memset(&ngs->st, 0, sizeof(ngs->st));
    ptmr_reset(&ngs->fwdtree_perf);
    ptmr_start(&ngs->fwdtree_perf);

    /* Reset backpointer table. */
    ngs->bpidx = 0;
    ngs->bss_head = 0;

    /* Reset word lattice. */
    for (i = 0; i < n_words; ++i)
        ngs->word_lat_idx[i] = NO_BP;

    /* Reset active HMM and word lists. */
    ngs->n_active_chan[0] = ngs->n_active_chan[1] = 0;
    ngs->n_active_word[0] = ngs->n_active_word[1] = 0;

    /* Reset scores. */
    ngs->best_score = 0;
    ngs->renormalized = 0;

    /* Reset other stuff. */
    for (i = 0; i < n_words; i++)
        ngs->last_ltrans[i].sf = -1;
    ngs->n_frame = 0;

    /* Clear the hypothesis string. */
    ckd_free(base->hyp_str);
    base->hyp_str = NULL;

    /* Reset the permanently allocated single-phone words, since they
     * may have junk left over in them from FWDFLAT. */
    for (i = 0; i < ngs->n_1ph_words; i++) {
        w = ngs->single_phone_wid[i];
        rhmm = (root_chan_t *) ngs->word_chan[w];
        hmm_clear(&rhmm->hmm);
    }

    /* Start search with <s>; word_chan[<s>] is permanently allocated */
    rhmm = (root_chan_t *) ngs->word_chan[dict_startwid(ps_search_dict(ngs))];
    hmm_clear(&rhmm->hmm);
    hmm_enter(&rhmm->hmm, 0, NO_BP, 0);
}

/*
 * Mark the active senones for all senones belonging to channels that are active in the
 * current frame.
 */
static void
compute_sen_active(ngram_search_t *ngs, int frame_idx)
{
    root_chan_t *rhmm;
    chan_t *hmm, **acl;
    int32 i, w, *awl;

    acmod_clear_active(ps_search_acmod(ngs));

    /* Flag active senones for root channels */
    for (i = ngs->n_root_chan, rhmm = ngs->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx)
            acmod_activate_hmm(ps_search_acmod(ngs), &rhmm->hmm);
    }

    /* Flag active senones for nonroot channels in HMM tree */
    i = ngs->n_active_chan[frame_idx & 0x1];
    acl = ngs->active_chan_list[frame_idx & 0x1];
    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++)) {
        acmod_activate_hmm(ps_search_acmod(ngs), &hmm->hmm);
    }

    /* Flag active senones for individual word channels */
    i = ngs->n_active_word[frame_idx & 0x1];
    awl = ngs->active_word_list[frame_idx & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        for (hmm = ngs->word_chan[w]; hmm; hmm = hmm->next) {
            acmod_activate_hmm(ps_search_acmod(ngs), &hmm->hmm);
        }
    }
    for (i = 0; i < ngs->n_1ph_words; i++) {
        w = ngs->single_phone_wid[i];
        rhmm = (root_chan_t *) ngs->word_chan[w];

        if (hmm_frame(&rhmm->hmm) == frame_idx)
            acmod_activate_hmm(ps_search_acmod(ngs), &rhmm->hmm);
    }
}

static void
renormalize_scores(ngram_search_t *ngs, int frame_idx, int32 norm)
{
    root_chan_t *rhmm;
    chan_t *hmm, **acl;
    int32 i, w, *awl;

    /* Renormalize root channels */
    for (i = ngs->n_root_chan, rhmm = ngs->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_normalize(&rhmm->hmm, norm);
        }
    }

    /* Renormalize nonroot channels in HMM tree */
    i = ngs->n_active_chan[frame_idx & 0x1];
    acl = ngs->active_chan_list[frame_idx & 0x1];
    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++)) {
        hmm_normalize(&hmm->hmm, norm);
    }

    /* Renormalize individual word channels */
    i = ngs->n_active_word[frame_idx & 0x1];
    awl = ngs->active_word_list[frame_idx & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        for (hmm = ngs->word_chan[w]; hmm; hmm = hmm->next) {
            hmm_normalize(&hmm->hmm, norm);
        }
    }
    for (i = 0; i < ngs->n_1ph_words; i++) {
        w = ngs->single_phone_wid[i];
        rhmm = (root_chan_t *) ngs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_normalize(&rhmm->hmm, norm);
        }
    }

    ngs->renormalized = TRUE;
}

static int32
eval_root_chan(ngram_search_t *ngs, int frame_idx)
{
    root_chan_t *rhmm;
    int32 i, bestscore;

    bestscore = WORST_SCORE;
    for (i = ngs->n_root_chan, rhmm = ngs->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            int32 score = chan_v_eval(rhmm);
            if (score BETTER_THAN bestscore)
                bestscore = score;
            ++ngs->st.n_root_chan_eval;
        }
    }
    return (bestscore);
}

static int32
eval_nonroot_chan(ngram_search_t *ngs, int frame_idx)
{
    chan_t *hmm, **acl;
    int32 i, bestscore;

    i = ngs->n_active_chan[frame_idx & 0x1];
    acl = ngs->active_chan_list[frame_idx & 0x1];
    bestscore = WORST_SCORE;
    ngs->st.n_nonroot_chan_eval += i;

    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++)) {
        int32 score = chan_v_eval(hmm);
        assert(hmm_frame(&hmm->hmm) == frame_idx);
        if (score BETTER_THAN bestscore)
            bestscore = score;
    }

    return bestscore;
}

static int32
eval_word_chan(ngram_search_t *ngs, int frame_idx)
{
    root_chan_t *rhmm;
    chan_t *hmm;
    int32 i, w, bestscore, *awl, j, k;

    k = 0;
    bestscore = WORST_SCORE;
    awl = ngs->active_word_list[frame_idx & 0x1];

    i = ngs->n_active_word[frame_idx & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        assert(bitvec_is_set(ngs->word_active, w));
        bitvec_clear(ngs->word_active, w);
        assert(ngs->word_chan[w] != NULL);

        for (hmm = ngs->word_chan[w]; hmm; hmm = hmm->next) {
            int32 score;

            assert(hmm_frame(&hmm->hmm) == frame_idx);
            score = chan_v_eval(hmm);
            /*printf("eval word chan %d score %d\n", w, score); */

            if (score BETTER_THAN bestscore)
                bestscore = score;

            k++;
        }
    }

    /* Similarly for statically allocated single-phone words */
    j = 0;
    for (i = 0; i < ngs->n_1ph_words; i++) {
        int32 score;

        w = ngs->single_phone_wid[i];
        rhmm = (root_chan_t *) ngs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) < frame_idx)
            continue;

        score = chan_v_eval(rhmm);
        /* printf("eval 1ph word chan %d score %d\n", w, score); */
        if (score BETTER_THAN bestscore && w != ps_search_finish_wid(ngs))
            bestscore = score;

        j++;
    }

    ngs->st.n_last_chan_eval += k + j;
    ngs->st.n_nonroot_chan_eval += k + j;
    ngs->st.n_word_lastchan_eval +=
        ngs->n_active_word[frame_idx & 0x1] + j;

    return bestscore;
}

static int32
evaluate_channels(ngram_search_t *ngs, int16 const *senone_scores, int frame_idx)
{
    int32 bs;

    hmm_context_set_senscore(ngs->hmmctx, senone_scores);
    ngs->best_score = eval_root_chan(ngs, frame_idx);
    if ((bs = eval_nonroot_chan(ngs, frame_idx)) BETTER_THAN ngs->best_score)
        ngs->best_score = bs;
    if ((bs = eval_word_chan(ngs, frame_idx)) BETTER_THAN ngs->best_score)
        ngs->best_score = bs;
    ngs->last_phone_best_score = bs;

    return ngs->best_score;
}

/*
 * Prune currently active root channels for next frame.  Also, perform exit
 * transitions out of them and activate successors.
 * score[] of pruned root chan set to WORST_SCORE elsewhere.
 */
static void
prune_root_chan(ngram_search_t *ngs, int frame_idx)
{
    root_chan_t *rhmm;
    chan_t *hmm;
    int32 i, nf, w;
    int32 thresh, newphone_thresh, lastphn_thresh, newphone_score;
    chan_t **nacl;              /* next active list */
    lastphn_cand_t *candp;
    phone_loop_search_t *pls;

    nf = frame_idx + 1;
    thresh = ngs->best_score + ngs->dynamic_beam;
    newphone_thresh = ngs->best_score + ngs->pbeam;
    lastphn_thresh = ngs->best_score + ngs->lpbeam;
    nacl = ngs->active_chan_list[nf & 0x1];
    pls = (phone_loop_search_t *)ps_search_lookahead(ngs);

    for (i = 0, rhmm = ngs->root_chan; i < ngs->n_root_chan; i++, rhmm++) {
        E_DEBUG(3,("Root channel %d frame %d score %d thresh %d\n",
                   i, hmm_frame(&rhmm->hmm), hmm_bestscore(&rhmm->hmm), thresh));
        /* First check if this channel was active in current frame */
        if (hmm_frame(&rhmm->hmm) < frame_idx)
            continue;

        if (hmm_bestscore(&rhmm->hmm) BETTER_THAN thresh) {
            hmm_frame(&rhmm->hmm) = nf;  /* rhmm will be active in next frame */
            E_DEBUG(3,("Preserving root channel %d score %d\n", i, hmm_bestscore(&rhmm->hmm)));
            /* transitions out of this root channel */
            /* transition to all next-level channels in the HMM tree */
            newphone_score = hmm_out_score(&rhmm->hmm) + ngs->pip;
            if (pls != NULL || newphone_score BETTER_THAN newphone_thresh) {
                for (hmm = rhmm->next; hmm; hmm = hmm->alt) {
                    int32 pl_newphone_score = newphone_score
                        + phone_loop_search_score(pls, hmm->ciphone);
                    if (pl_newphone_score BETTER_THAN newphone_thresh) {
                        if ((hmm_frame(&hmm->hmm) < frame_idx)
                            || (pl_newphone_score BETTER_THAN hmm_in_score(&hmm->hmm))) {
                            hmm_enter(&hmm->hmm, pl_newphone_score,
                                      hmm_out_history(&rhmm->hmm), nf);
                            *(nacl++) = hmm;
                        }
                    }
                }
            }

            /*
             * Transition to last phone of all words for which this is the
             * penultimate phone (the last phones may need multiple right contexts).
             * Remember to remove the temporary newword_penalty.
             */
            if (pls != NULL || newphone_score BETTER_THAN lastphn_thresh) {
                for (w = rhmm->penult_phn_wid; w >= 0;
                     w = ngs->homophone_set[w]) {
                    int32 pl_newphone_score = newphone_score
                        + phone_loop_search_score
                        (pls, dict_last_phone(ps_search_dict(ngs),w));
                    E_DEBUG(3,("word %s newphone_score %d\n", dict_wordstr(ps_search_dict(ngs), w), newphone_score));
                    if (pl_newphone_score BETTER_THAN lastphn_thresh) {
                        candp = ngs->lastphn_cand + ngs->n_lastphn_cand;
                        ngs->n_lastphn_cand++;
                        candp->wid = w;
                        candp->score =
                            pl_newphone_score - ngs->nwpen;
                        candp->bp = hmm_out_history(&rhmm->hmm);
                    }
                }
            }
        }
    }
    ngs->n_active_chan[nf & 0x1] = nacl - ngs->active_chan_list[nf & 0x1];
}

/*
 * Prune currently active nonroot channels in HMM tree for next frame.  Also, perform
 * exit transitions out of such channels and activate successors.
 */
static void
prune_nonroot_chan(ngram_search_t *ngs, int frame_idx)
{
    chan_t *hmm, *nexthmm;
    int32 nf, w, i;
    int32 thresh, newphone_thresh, lastphn_thresh, newphone_score;
    chan_t **acl, **nacl;       /* active list, next active list */
    lastphn_cand_t *candp;
    phone_loop_search_t *pls;

    nf = frame_idx + 1;

    thresh = ngs->best_score + ngs->dynamic_beam;
    newphone_thresh = ngs->best_score + ngs->pbeam;
    lastphn_thresh = ngs->best_score + ngs->lpbeam;
    pls = (phone_loop_search_t *)ps_search_lookahead(ngs);

    acl = ngs->active_chan_list[frame_idx & 0x1];   /* currently active HMMs in tree */
    nacl = ngs->active_chan_list[nf & 0x1] + ngs->n_active_chan[nf & 0x1];

    for (i = ngs->n_active_chan[frame_idx & 0x1], hmm = *(acl++); i > 0;
         --i, hmm = *(acl++)) {
        assert(hmm_frame(&hmm->hmm) >= frame_idx);

        if (hmm_bestscore(&hmm->hmm) BETTER_THAN thresh) {
            /* retain this channel in next frame */
            if (hmm_frame(&hmm->hmm) != nf) {
                hmm_frame(&hmm->hmm) = nf;
                *(nacl++) = hmm;
            }

            /* transition to all next-level channel in the HMM tree */
            newphone_score = hmm_out_score(&hmm->hmm) + ngs->pip;
            if (pls != NULL || newphone_score BETTER_THAN newphone_thresh) {
                for (nexthmm = hmm->next; nexthmm; nexthmm = nexthmm->alt) {
                    int32 pl_newphone_score = newphone_score
                        + phone_loop_search_score(pls, nexthmm->ciphone);
                    if ((pl_newphone_score BETTER_THAN newphone_thresh)
                        && ((hmm_frame(&nexthmm->hmm) < frame_idx)
                            || (pl_newphone_score
                                BETTER_THAN hmm_in_score(&nexthmm->hmm)))) {
                        if (hmm_frame(&nexthmm->hmm) != nf) {
                            /* Keep this HMM on the active list */
                            *(nacl++) = nexthmm;
                        }
                        hmm_enter(&nexthmm->hmm, pl_newphone_score,
                                  hmm_out_history(&hmm->hmm), nf);
                    }
                }
            }

            /*
             * Transition to last phone of all words for which this is the
             * penultimate phone (the last phones may need multiple right contexts).
             * Remember to remove the temporary newword_penalty.
             */
            if (pls != NULL || newphone_score BETTER_THAN lastphn_thresh) {
                for (w = hmm->info.penult_phn_wid; w >= 0;
                     w = ngs->homophone_set[w]) {
                    int32 pl_newphone_score = newphone_score
                        + phone_loop_search_score
                        (pls, dict_last_phone(ps_search_dict(ngs),w));
                    if (pl_newphone_score BETTER_THAN lastphn_thresh) {
                        candp = ngs->lastphn_cand + ngs->n_lastphn_cand;
                        ngs->n_lastphn_cand++;
                        candp->wid = w;
                        candp->score =
                            pl_newphone_score - ngs->nwpen;
                        candp->bp = hmm_out_history(&hmm->hmm);
                    }
                }
            }
        }
        else if (hmm_frame(&hmm->hmm) != nf) {
            hmm_clear(&hmm->hmm);
        }
    }
    ngs->n_active_chan[nf & 0x1] = nacl - ngs->active_chan_list[nf & 0x1];
}

/*
 * Execute the transition into the last phone for all candidates words emerging from
 * the HMM tree.  Attach LM scores to such transitions.
 * (Executed after pruning root and non-root, but before pruning word-chan.)
 */
static void
last_phone_transition(ngram_search_t *ngs, int frame_idx)
{
    int32 i, j, k, nf, bp, bpend, w;
    lastphn_cand_t *candp;
    int32 *nawl;
    int32 thresh;
    int32 bestscore, dscr;
    chan_t *hmm;
    bptbl_t *bpe;
    int32 n_cand_sf = 0;

    nf = frame_idx + 1;
    nawl = ngs->active_word_list[nf & 0x1];
    ngs->st.n_lastphn_cand_utt += ngs->n_lastphn_cand;

    /* For each candidate word (entering its last phone) */
    /* If best LM score and bp for candidate known use it, else sort cands by startfrm */
    for (i = 0, candp = ngs->lastphn_cand; i < ngs->n_lastphn_cand; i++, candp++) {
        int32 start_score;

        /* This can happen if recognition fails. */
        if (candp->bp == -1)
            continue;
        /* Backpointer entry for it. */
        bpe = &(ngs->bp_table[candp->bp]);

        /* Subtract starting score for candidate, leave it with only word score */
        start_score = ngram_search_exit_score
            (ngs, bpe, dict_first_phone(ps_search_dict(ngs), candp->wid));
        assert(start_score BETTER_THAN WORST_SCORE);
        candp->score -= start_score;

        /*
         * If this candidate not occurred in an earlier frame, prepare for finding
         * best transition score into last phone; sort by start frame.
         */
        /* i.e. if we don't have an entry in last_ltrans for this
         * <word,sf>, then create one */
        if (ngs->last_ltrans[candp->wid].sf != bpe->frame + 1) {
            /* Look for an entry in cand_sf matching the backpointer
             * for this candidate. */
            for (j = 0; j < n_cand_sf; j++) {
                if (ngs->cand_sf[j].bp_ef == bpe->frame)
                    break;
            }
            /* Oh, we found one, so chain onto it. */
            if (j < n_cand_sf)
                candp->next = ngs->cand_sf[j].cand;
            else {
                /* Nope, let's make a new one, allocating cand_sf if necessary. */
                if (n_cand_sf >= ngs->cand_sf_alloc) {
                    if (ngs->cand_sf_alloc == 0) {
                        ngs->cand_sf =
                            ckd_calloc(CAND_SF_ALLOCSIZE,
                                       sizeof(*ngs->cand_sf));
                        ngs->cand_sf_alloc = CAND_SF_ALLOCSIZE;
                    }
                    else {
                        ngs->cand_sf_alloc += CAND_SF_ALLOCSIZE;
                        ngs->cand_sf = ckd_realloc(ngs->cand_sf,
                                                   ngs->cand_sf_alloc
                                                   * sizeof(*ngs->cand_sf));
                        E_INFO("cand_sf[] increased to %d entries\n",
                               ngs->cand_sf_alloc);
                    }
                }

                /* Use the newly created cand_sf. */
                j = n_cand_sf++;
                candp->next = -1; /* End of the chain. */
                ngs->cand_sf[j].bp_ef = bpe->frame;
            }
            /* Update it to point to this candidate. */
            ngs->cand_sf[j].cand = i;

            ngs->last_ltrans[candp->wid].dscr = WORST_SCORE;
            ngs->last_ltrans[candp->wid].sf = bpe->frame + 1;
        }
    }

    /* Compute best LM score and bp for new cands entered in the sorted lists above */
    for (i = 0; i < n_cand_sf; i++) {
        /* For the i-th unique end frame... */
        bp = ngs->bp_table_idx[ngs->cand_sf[i].bp_ef];
        bpend = ngs->bp_table_idx[ngs->cand_sf[i].bp_ef + 1];
        for (bpe = &(ngs->bp_table[bp]); bp < bpend; bp++, bpe++) {
            if (!bpe->valid)
                continue;
            /* For each candidate at the start frame find bp->cand transition-score */
            for (j = ngs->cand_sf[i].cand; j >= 0; j = candp->next) {
                int32 n_used;
                candp = &(ngs->lastphn_cand[j]);
                dscr = 
                    ngram_search_exit_score
                    (ngs, bpe, dict_first_phone(ps_search_dict(ngs), candp->wid));
                if (dscr BETTER_THAN WORST_SCORE) {
                    assert(!dict_filler_word(ps_search_dict(ngs), candp->wid));
                    dscr += ngram_tg_score(ngs->lmset,
                                           dict_basewid(ps_search_dict(ngs), candp->wid),
                                           bpe->real_wid,
                                           bpe->prev_real_wid,
                                           &n_used)>>SENSCR_SHIFT;
                }

                if (dscr BETTER_THAN ngs->last_ltrans[candp->wid].dscr) {
                    ngs->last_ltrans[candp->wid].dscr = dscr;
                    ngs->last_ltrans[candp->wid].bp = bp;
                }
            }
        }
    }

    /* Update best transitions for all candidates; also update best lastphone score */
    bestscore = ngs->last_phone_best_score;
    for (i = 0, candp = ngs->lastphn_cand; i < ngs->n_lastphn_cand; i++, candp++) {
        candp->score += ngs->last_ltrans[candp->wid].dscr;
        candp->bp = ngs->last_ltrans[candp->wid].bp;

        if (candp->score BETTER_THAN bestscore)
            bestscore = candp->score;
    }
    ngs->last_phone_best_score = bestscore;

    /* At this pt, we know the best entry score (with LM component) for all candidates */
    thresh = bestscore + ngs->lponlybeam;
    for (i = ngs->n_lastphn_cand, candp = ngs->lastphn_cand; i > 0; --i, candp++) {
        if (candp->score BETTER_THAN thresh) {
            w = candp->wid;

            ngram_search_alloc_all_rc(ngs, w);

            k = 0;
            for (hmm = ngs->word_chan[w]; hmm; hmm = hmm->next) {
                if ((hmm_frame(&hmm->hmm) < frame_idx)
                    || (candp->score BETTER_THAN hmm_in_score(&hmm->hmm))) {
                    assert(hmm_frame(&hmm->hmm) != nf);
                    hmm_enter(&hmm->hmm,
                              candp->score, candp->bp, nf);
                    k++;
                }
            }
            if (k > 0) {
                assert(bitvec_is_clear(ngs->word_active, w));
                assert(!dict_is_single_phone(ps_search_dict(ngs), w));
                *(nawl++) = w;
                bitvec_set(ngs->word_active, w);
            }
        }
    }
    ngs->n_active_word[nf & 0x1] = nawl - ngs->active_word_list[nf & 0x1];
}

/*
 * Prune currently active word channels for next frame.  Also, perform exit
 * transitions out of such channels and active successors.
 */
static void
prune_word_chan(ngram_search_t *ngs, int frame_idx)
{
    root_chan_t *rhmm;
    chan_t *hmm, *thmm;
    chan_t **phmmp;             /* previous HMM-pointer */
    int32 nf, w, i, k;
    int32 newword_thresh, lastphn_thresh;
    int32 *awl, *nawl;

    nf = frame_idx + 1;
    newword_thresh = ngs->last_phone_best_score + ngs->wbeam;
    lastphn_thresh = ngs->last_phone_best_score + ngs->lponlybeam;

    awl = ngs->active_word_list[frame_idx & 0x1];
    nawl = ngs->active_word_list[nf & 0x1] + ngs->n_active_word[nf & 0x1];

    /* Dynamically allocated last channels of multi-phone words */
    for (i = ngs->n_active_word[frame_idx & 0x1], w = *(awl++); i > 0;
         --i, w = *(awl++)) {
        k = 0;
        phmmp = &(ngs->word_chan[w]);
        for (hmm = ngs->word_chan[w]; hmm; hmm = thmm) {
            assert(hmm_frame(&hmm->hmm) >= frame_idx);

            thmm = hmm->next;
            if (hmm_bestscore(&hmm->hmm) BETTER_THAN lastphn_thresh) {
                /* retain this channel in next frame */
                hmm_frame(&hmm->hmm) = nf;
                k++;
                phmmp = &(hmm->next);

                /* Could if ((! skip_alt_frm) || (frame_idx & 0x1)) the following */
                if (hmm_out_score(&hmm->hmm) BETTER_THAN newword_thresh) {
                    /* can exit channel and recognize word */
                    ngram_search_save_bp(ngs, frame_idx, w,
                                 hmm_out_score(&hmm->hmm),
                                 hmm_out_history(&hmm->hmm),
                                 hmm->info.rc_id);
                }
            }
            else if (hmm_frame(&hmm->hmm) == nf) {
                phmmp = &(hmm->next);
            }
            else {
                hmm_deinit(&hmm->hmm);
                listelem_free(ngs->chan_alloc, hmm);
                *phmmp = thmm;
            }
        }
        if ((k > 0) && (bitvec_is_clear(ngs->word_active, w))) {
            assert(!dict_is_single_phone(ps_search_dict(ngs), w));
            *(nawl++) = w;
            bitvec_set(ngs->word_active, w);
        }
    }
    ngs->n_active_word[nf & 0x1] = nawl - ngs->active_word_list[nf & 0x1];

    /*
     * Prune permanently allocated single-phone channels.
     * NOTES: score[] of pruned channels set to WORST_SCORE elsewhere.
     */
    for (i = 0; i < ngs->n_1ph_words; i++) {
        w = ngs->single_phone_wid[i];
        rhmm = (root_chan_t *) ngs->word_chan[w];
        E_DEBUG(3,("Single phone word %s frame %d score %d thresh %d outscore %d nwthresh %d\n",
                   dict_wordstr(ps_search_dict(ngs),w),
                   hmm_frame(&rhmm->hmm), hmm_bestscore(&rhmm->hmm),
                   lastphn_thresh, hmm_out_score(&rhmm->hmm), newword_thresh));
        if (hmm_frame(&rhmm->hmm) < frame_idx)
            continue;
        if (hmm_bestscore(&rhmm->hmm) BETTER_THAN lastphn_thresh) {
            hmm_frame(&rhmm->hmm) = nf;

            /* Could if ((! skip_alt_frm) || (frame_idx & 0x1)) the following */
            if (hmm_out_score(&rhmm->hmm) BETTER_THAN newword_thresh) {
                E_DEBUG(4,("Exiting single phone word %s with %d > %d, %d\n",
                           dict_wordstr(ps_search_dict(ngs),w),
                           hmm_out_score(&rhmm->hmm),
                           lastphn_thresh, newword_thresh));
                ngram_search_save_bp(ngs, frame_idx, w,
                             hmm_out_score(&rhmm->hmm),
                             hmm_out_history(&rhmm->hmm), 0);
            }
        }
    }
}

static void
prune_channels(ngram_search_t *ngs, int frame_idx)
{
    /* Clear last phone candidate list. */
    ngs->n_lastphn_cand = 0;
    /* Set the dynamic beam based on maxhmmpf here. */
    ngs->dynamic_beam = ngs->beam;
    if (ngs->maxhmmpf != -1
        && ngs->st.n_root_chan_eval + ngs->st.n_nonroot_chan_eval > ngs->maxhmmpf) {
        /* Build a histogram to approximately prune them. */
        int32 bins[256], bw, nhmms, i;
        root_chan_t *rhmm;
        chan_t **acl, *hmm;

        /* Bins go from zero (best score) to edge of beam. */
        bw = -ngs->beam / 256;
        memset(bins, 0, sizeof(bins));
        /* For each active root channel. */
        for (i = 0, rhmm = ngs->root_chan; i < ngs->n_root_chan; i++, rhmm++) {
            int32 b;

            /* Put it in a bin according to its bestscore. */
            b = (ngs->best_score - hmm_bestscore(&rhmm->hmm)) / bw;
            if (b >= 256)
                b = 255;
            ++bins[b];
        }
        /* For each active non-root channel. */
        acl = ngs->active_chan_list[frame_idx & 0x1];       /* currently active HMMs in tree */
        for (i = ngs->n_active_chan[frame_idx & 0x1], hmm = *(acl++);
             i > 0; --i, hmm = *(acl++)) {
            int32 b;

            /* Put it in a bin according to its bestscore. */
            b = (ngs->best_score - hmm_bestscore(&hmm->hmm)) / bw;
            if (b >= 256)
                b = 255;
            ++bins[b];
        }
        /* Walk down the bins to find the new beam. */
        for (i = nhmms = 0; i < 256; ++i) {
            nhmms += bins[i];
            if (nhmms > ngs->maxhmmpf)
                break;
        }
        ngs->dynamic_beam = -(i * bw);
    }

    prune_root_chan(ngs, frame_idx);
    prune_nonroot_chan(ngs, frame_idx);
    last_phone_transition(ngs, frame_idx);
    prune_word_chan(ngs, frame_idx);
}

/*
 * Limit the number of word exits in each frame to maxwpf.  And also limit the number of filler
 * words to 1.
 */
static void
bptable_maxwpf(ngram_search_t *ngs, int frame_idx)
{
    int32 bp, n;
    int32 bestscr, worstscr;
    bptbl_t *bpe, *bestbpe, *worstbpe;

    /* Don't prune if no pruing. */
    if (ngs->maxwpf == -1 || ngs->maxwpf == ps_search_n_words(ngs))
        return;

    /* Allow only one filler word exit (the best) per frame */
    bestscr = (int32) 0x80000000;
    bestbpe = NULL;
    n = 0;
    for (bp = ngs->bp_table_idx[frame_idx]; bp < ngs->bpidx; bp++) {
        bpe = &(ngs->bp_table[bp]);
        if (dict_filler_word(ps_search_dict(ngs), bpe->wid)) {
            if (bpe->score BETTER_THAN bestscr) {
                bestscr = bpe->score;
                bestbpe = bpe;
            }
            bpe->valid = FALSE;
            n++;                /* No. of filler words */
        }
    }
    /* Restore bestbpe to valid state */
    if (bestbpe != NULL) {
        bestbpe->valid = TRUE;
        --n;
    }

    /* Allow up to maxwpf best entries to survive; mark the remaining with valid = 0 */
    n = (ngs->bpidx
         - ngs->bp_table_idx[frame_idx]) - n;  /* No. of entries after limiting fillers */
    for (; n > ngs->maxwpf; --n) {
        /* Find worst BPTable entry */
        worstscr = (int32) 0x7fffffff;
        worstbpe = NULL;
        for (bp = ngs->bp_table_idx[frame_idx]; (bp < ngs->bpidx); bp++) {
            bpe = &(ngs->bp_table[bp]);
            if (bpe->valid && (bpe->score WORSE_THAN worstscr)) {
                worstscr = bpe->score;
                worstbpe = bpe;
            }
        }
        /* FIXME: Don't panic! */
        if (worstbpe == NULL)
            E_FATAL("PANIC: No worst BPtable entry remaining\n");
        worstbpe->valid = FALSE;
    }
}

static void
word_transition(ngram_search_t *ngs, int frame_idx)
{
    int32 i, k, bp, w, nf;
    int32 rc;
    int32 thresh, newscore;
    bptbl_t *bpe;
    root_chan_t *rhmm;
    struct bestbp_rc_s *bestbp_rc_ptr;
    phone_loop_search_t *pls;
    dict_t *dict = ps_search_dict(ngs);
    dict2pid_t *d2p = ps_search_dict2pid(ngs);

    /*
     * Transition to start of new word instances (HMM tree roots); but only if words
     * other than </s> finished here.
     * But, first, find the best starting score for each possible right context phone.
     */
    for (i = bin_mdef_n_ciphone(ps_search_acmod(ngs)->mdef) - 1; i >= 0; --i)
        ngs->bestbp_rc[i].score = WORST_SCORE;
    k = 0;
    pls = (phone_loop_search_t *)ps_search_lookahead(ngs);
    /* Ugh, this is complicated.  Scan all word exits for this frame
     * (they have already been created by prune_word_chan()). */
    for (bp = ngs->bp_table_idx[frame_idx]; bp < ngs->bpidx; bp++) {
        bpe = &(ngs->bp_table[bp]);
        ngs->word_lat_idx[bpe->wid] = NO_BP;

        if (bpe->wid == ps_search_finish_wid(ngs))
            continue;
        k++;

        /* DICT2PID */
        /* Array of HMM scores corresponding to all the possible right
         * context expansions of the final phone.  It's likely that a
         * lot of these are going to be missing, actually. */
        if (bpe->last2_phone == -1) { /* implies s_idx == -1 */
            /* No right context expansion. */
            for (rc = 0; rc < bin_mdef_n_ciphone(ps_search_acmod(ngs)->mdef); ++rc) {
                if (bpe->score BETTER_THAN ngs->bestbp_rc[rc].score) {
                    E_DEBUG(4,("bestbp_rc[0] = %d lc %d\n",
                               bpe->score, bpe->last_phone));
                    ngs->bestbp_rc[rc].score = bpe->score;
                    ngs->bestbp_rc[rc].path = bp;
                    ngs->bestbp_rc[rc].lc = bpe->last_phone;
                }
            }
        }
        else {
            xwdssid_t *rssid = dict2pid_rssid(d2p, bpe->last_phone, bpe->last2_phone);
            int32 *rcss = &(ngs->bscore_stack[bpe->s_idx]);
            for (rc = 0; rc < bin_mdef_n_ciphone(ps_search_acmod(ngs)->mdef); ++rc) {
                if (rcss[rssid->cimap[rc]] BETTER_THAN ngs->bestbp_rc[rc].score) {
                    E_DEBUG(4,("bestbp_rc[%d] = %d lc %d\n",
                               rc, rcss[rssid->cimap[rc]], bpe->last_phone));
                    ngs->bestbp_rc[rc].score = rcss[rssid->cimap[rc]];
                    ngs->bestbp_rc[rc].path = bp;
                    ngs->bestbp_rc[rc].lc = bpe->last_phone;
                }
            }
        }
    }
    if (k == 0)
        return;

    nf = frame_idx + 1;
    thresh = ngs->best_score + ngs->dynamic_beam;
    /*
     * Hypothesize successors to words finished in this frame.
     * Main dictionary, multi-phone words transition to HMM-trees roots.
     */
    for (i = ngs->n_root_chan, rhmm = ngs->root_chan; i > 0; --i, rhmm++) {
        bestbp_rc_ptr = &(ngs->bestbp_rc[rhmm->ciphone]);

        newscore = bestbp_rc_ptr->score + ngs->nwpen + ngs->pip
            + phone_loop_search_score(pls, rhmm->ciphone);
        if (newscore BETTER_THAN thresh) {
            if ((hmm_frame(&rhmm->hmm) < frame_idx)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm, newscore,
                          bestbp_rc_ptr->path, nf);
                /* DICT2PID: Another place where mpx ssids are entered. */
                /* Look up the ssid to use when entering this mpx triphone. */
                hmm_mpx_ssid(&rhmm->hmm, 0) =
                    dict2pid_ldiph_lc(d2p, rhmm->ciphone, rhmm->ci2phone, bestbp_rc_ptr->lc);
                assert(hmm_mpx_ssid(&rhmm->hmm, 0) != BAD_SSID);
            }
        }
    }

    /*
     * Single phone words; no right context for these.  Cannot use bestbp_rc as
     * LM scores have to be included.  First find best transition to these words.
     */
    for (i = 0; i < ngs->n_1ph_LMwords; i++) {
        w = ngs->single_phone_wid[i];
        ngs->last_ltrans[w].dscr = (int32) 0x80000000;
    }
    for (bp = ngs->bp_table_idx[frame_idx]; bp < ngs->bpidx; bp++) {
        bpe = &(ngs->bp_table[bp]);
        if (!bpe->valid)
            continue;

        for (i = 0; i < ngs->n_1ph_LMwords; i++) {
            int32 n_used;
            w = ngs->single_phone_wid[i];
            newscore = ngram_search_exit_score
                (ngs, bpe, dict_first_phone(dict, w));
            E_DEBUG(4, ("initial newscore for %s: %d\n",
                        dict_wordstr(dict, w), newscore));
            if (newscore != WORST_SCORE)
                newscore += ngram_tg_score(ngs->lmset,
                                           dict_basewid(dict, w),
                                           bpe->real_wid,
                                           bpe->prev_real_wid,
                                           &n_used)>>SENSCR_SHIFT;

            /* FIXME: Not sure how WORST_SCORE could be better, but it
             * apparently happens. */
            if (newscore BETTER_THAN ngs->last_ltrans[w].dscr) {
                ngs->last_ltrans[w].dscr = newscore;
                ngs->last_ltrans[w].bp = bp;
            }
        }
    }

    /* Now transition to in-LM single phone words */
    for (i = 0; i < ngs->n_1ph_LMwords; i++) {
        w = ngs->single_phone_wid[i];
        /* Never transition into the start word (for one thing, it is
           a non-event in the language model.) */
        if (w == dict_startwid(ps_search_dict(ngs)))
            continue;
        rhmm = (root_chan_t *) ngs->word_chan[w];
        newscore = ngs->last_ltrans[w].dscr + ngs->pip
            + phone_loop_search_score(pls, rhmm->ciphone);
        if (newscore BETTER_THAN thresh) {
            bpe = ngs->bp_table + ngs->last_ltrans[w].bp;
            if ((hmm_frame(&rhmm->hmm) < frame_idx)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm,
                          newscore, ngs->last_ltrans[w].bp, nf);
                /* DICT2PID: another place where mpx ssids are entered. */
                /* Look up the ssid to use when entering this mpx triphone. */
                hmm_mpx_ssid(&rhmm->hmm, 0) =
                    dict2pid_ldiph_lc(d2p, rhmm->ciphone, rhmm->ci2phone,
                                      dict_last_phone(dict, bpe->wid));
                assert(hmm_mpx_ssid(&rhmm->hmm, 0) != BAD_SSID);
            }
        }
    }

    /* Remaining words: <sil>, noise words.  No mpx for these! */
    w = ps_search_silence_wid(ngs);
    rhmm = (root_chan_t *) ngs->word_chan[w];
    bestbp_rc_ptr = &(ngs->bestbp_rc[ps_search_acmod(ngs)->mdef->sil]);
    newscore = bestbp_rc_ptr->score + ngs->silpen + ngs->pip
        + phone_loop_search_score(pls, rhmm->ciphone);
    if (newscore BETTER_THAN thresh) {
        if ((hmm_frame(&rhmm->hmm) < frame_idx)
            || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
            hmm_enter(&rhmm->hmm,
                      newscore, bestbp_rc_ptr->path, nf);
        }
    }
    for (w = dict_filler_start(dict); w <= dict_filler_end(dict); w++) {
        if (w == ps_search_silence_wid(ngs))
            continue;
        /* Never transition into the start word (for one thing, it is
           a non-event in the language model.) */
        if (w == dict_startwid(ps_search_dict(ngs)))
            continue;
        rhmm = (root_chan_t *) ngs->word_chan[w];
        /* If this was not actually a single-phone word, rhmm will be NULL. */
        if (rhmm == NULL)
            continue;
        newscore = bestbp_rc_ptr->score + ngs->fillpen + ngs->pip
            + phone_loop_search_score(pls, rhmm->ciphone);
        if (newscore BETTER_THAN thresh) {
            if ((hmm_frame(&rhmm->hmm) < frame_idx)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm,
                          newscore, bestbp_rc_ptr->path, nf);
            }
        }
    }
}

static void
deactivate_channels(ngram_search_t *ngs, int frame_idx)
{
    root_chan_t *rhmm;
    int i;

    /* Clear score[] of pruned root channels */
    for (i = ngs->n_root_chan, rhmm = ngs->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_clear(&rhmm->hmm);
        }
    }
    /* Clear score[] of pruned single-phone channels */
    for (i = 0; i < ngs->n_1ph_words; i++) {
        int32 w = ngs->single_phone_wid[i];
        rhmm = (root_chan_t *) ngs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_clear(&rhmm->hmm);
        }
    }
}

int
ngram_fwdtree_search(ngram_search_t *ngs, int frame_idx)
{
    int16 const *senscr;

    /* Activate our HMMs for the current frame if need be. */
    if (!ps_search_acmod(ngs)->compallsen)
        compute_sen_active(ngs, frame_idx);

    /* Compute GMM scores for the current frame. */
    if ((senscr = acmod_score(ps_search_acmod(ngs), &frame_idx)) == NULL)
        return 0;
    ngs->st.n_senone_active_utt += ps_search_acmod(ngs)->n_senone_active;

    /* Mark backpointer table for current frame. */
    ngram_search_mark_bptable(ngs, frame_idx);

    /* If the best score is equal to or worse than WORST_SCORE,
     * recognition has failed, don't bother to keep trying. */
    if (ngs->best_score == WORST_SCORE || ngs->best_score WORSE_THAN WORST_SCORE)
        return 0;
    /* Renormalize if necessary */
    if (ngs->best_score + (2 * ngs->beam) WORSE_THAN WORST_SCORE) {
        E_INFO("Renormalizing Scores at frame %d, best score %d\n",
               frame_idx, ngs->best_score);
        renormalize_scores(ngs, frame_idx, ngs->best_score);
    }

    /* Evaluate HMMs */
    evaluate_channels(ngs, senscr, frame_idx);
    /* Prune HMMs and do phone transitions. */
    prune_channels(ngs, frame_idx);
    /* Do absolute pruning on word exits. */
    bptable_maxwpf(ngs, frame_idx);
    /* Do word transitions. */
    word_transition(ngs, frame_idx);
    /* Deactivate pruned HMMs. */
    deactivate_channels(ngs, frame_idx);

    ++ngs->n_frame;
    /* Return the number of frames processed. */
    return 1;
}

void
ngram_fwdtree_finish(ngram_search_t *ngs)
{
    int32 i, w, cf, *awl;
    root_chan_t *rhmm;
    chan_t *hmm, **acl;

    /* This is the number of frames processed. */
    cf = ps_search_acmod(ngs)->output_frame;
    /* Add a mark in the backpointer table for one past the final frame. */
    ngram_search_mark_bptable(ngs, cf);

    /* Deactivate channels lined up for the next frame */
    /* First, root channels of HMM tree */
    for (i = ngs->n_root_chan, rhmm = ngs->root_chan; i > 0; --i, rhmm++) {
        hmm_clear(&rhmm->hmm);
    }

    /* nonroot channels of HMM tree */
    i = ngs->n_active_chan[cf & 0x1];
    acl = ngs->active_chan_list[cf & 0x1];
    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++)) {
        hmm_clear(&hmm->hmm);
    }

    /* word channels */
    i = ngs->n_active_word[cf & 0x1];
    awl = ngs->active_word_list[cf & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        /* Don't accidentally free single-phone words! */
        if (dict_is_single_phone(ps_search_dict(ngs), w))
            continue;
        bitvec_clear(ngs->word_active, w);
        if (ngs->word_chan[w] == NULL)
            continue;
        ngram_search_free_all_rc(ngs, w);
    }

    /*
     * The previous search code did a postprocessing of the
     * backpointer table here, but we will postpone this until it is
     * absolutely necessary, i.e. when generating a word graph.
     * Likewise we don't actually have to decide what the exit word is
     * until somebody requests a backtrace.
     */

    ptmr_stop(&ngs->fwdtree_perf);
    /* Print out some statistics. */
    if (cf > 0) {
        double n_speech = (double)(cf + 1)
            / cmd_ln_int32_r(ps_search_config(ngs), "-frate");
        E_INFO("%8d words recognized (%d/fr)\n",
               ngs->bpidx, (ngs->bpidx + (cf >> 1)) / (cf + 1));
        E_INFO("%8d senones evaluated (%d/fr)\n", ngs->st.n_senone_active_utt,
               (ngs->st.n_senone_active_utt + (cf >> 1)) / (cf + 1));
        E_INFO("%8d channels searched (%d/fr), %d 1st, %d last\n",
               ngs->st.n_root_chan_eval + ngs->st.n_nonroot_chan_eval,
               (ngs->st.n_root_chan_eval + ngs->st.n_nonroot_chan_eval) / (cf + 1),
               ngs->st.n_root_chan_eval, ngs->st.n_last_chan_eval);
        E_INFO("%8d words for which last channels evaluated (%d/fr)\n",
               ngs->st.n_word_lastchan_eval,
               ngs->st.n_word_lastchan_eval / (cf + 1));
        E_INFO("%8d candidate words for entering last phone (%d/fr)\n",
               ngs->st.n_lastphn_cand_utt, ngs->st.n_lastphn_cand_utt / (cf + 1));
        E_INFO("fwdtree %.2f CPU %.3f xRT\n",
               ngs->fwdtree_perf.t_cpu,
               ngs->fwdtree_perf.t_cpu / n_speech);
        E_INFO("fwdtree %.2f wall %.3f xRT\n",
               ngs->fwdtree_perf.t_elapsed,
               ngs->fwdtree_perf.t_elapsed / n_speech);
    }
    /* dump_bptable(ngs); */
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ngram_search_fwdtree.h Lexicon tree based Viterbi search.
 */

#ifndef __NGRAM_SEARCH_FWDTREE_H__
#define __NGRAM_SEARCH_FWDTREE_H__

/* SphinxBase headers. */

/* Local headers. */
#include "ngram_search.h"

/**
 * Initialize N-Gram search for fwdtree decoding.
 */
void ngram_fwdtree_init(ngram_search_t *ngs);

/**
 * Release memory associated with fwdtree decoding.
 */
void ngram_fwdtree_deinit(ngram_search_t *ngs);

/**
 * Rebuild search structures for updated language models.
 */
int ngram_fwdtree_reinit(ngram_search_t *ngs);

/**
 * Start fwdtree decoding for an utterance.
 */
void ngram_fwdtree_start(ngram_search_t *ngs);

/**
 * Search one frame forward in an utterance.
 *
 * @return Number of frames searched (either 0 or 1).
 */
int ngram_fwdtree_search(ngram_search_t *ngs, int frame_idx);

/**
 * Finish fwdtree decoding for an utterance.
 */
void ngram_fwdtree_finish(ngram_search_t *ngs);


#endif /* __NGRAM_SEARCH_FWDTREE_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ngram_search.h N-Gram based multi-pass search ("FBS")
 */

#ifndef __NGRAM_SEARCH_H__
#define __NGRAM_SEARCH_H__

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/ngram_model.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "pocketsphinx_internal.h"
#include "hmm.h"

/**
 * Lexical tree node data type.
 *
 * Not the first HMM for words, which multiplex HMMs based on
 * different left contexts.  This structure is used both in the
 * dynamic HMM tree structure and in the per-word last-phone right
 * context fanout.
 */
typedef struct chan_s {
    hmm_t hmm;                  /**< Basic HMM structure.  This *must* be first in
                                   the structure because chan_t and root_chan_t are
                                   sometimes used interchangeably */
    struct chan_s *next;	/**< first descendant of this channel; or, in the
				   case of the last phone of a word, the next
				   alternative right context channel */
    struct chan_s *alt;		/**< sibling; i.e., next descendant of parent HMM */

    int32    ciphone;		/**< ciphone for this node */
    union {
	int32 penult_phn_wid;	/**< list of words whose last phone follows this one;
				   this field indicates the first of the list; the
				   rest must be built up in a separate array.  Used
				   only within HMM tree.  -1 if none */
	int32 rc_id;		/**< right-context id for last phone of words */
    } info;
} chan_t;

/**
 * Lexical tree node data type for the first phone (root) of each dynamic HMM tree
 * structure.
 *
 * Each state may have a different parent static HMM.  Most fields are
 * similar to those in chan_t.
 */
typedef struct root_chan_s {
    hmm_t hmm;                  /**< Basic HMM structure.  This *must* be first in
                                   the structure because chan_t and root_chan_t are
                                   sometimes used interchangeably. */
    chan_t *next;		/**< first descendant of this channel */

    int32  penult_phn_wid;
    int32  this_phn_wid;	/**< list of words consisting of this single phone;
				   actually the first of the list, like penult_phn_wid;
				   -1 if none */
    int16    ciphone;		/**< first ciphone of this node; all words rooted at this
				   node begin with this ciphone */
    int16    ci2phone;		/**< second ciphone of this node; one root HMM for each
                                   unique right context */
} root_chan_t;

/**
 * Back pointer table (forward pass lattice; actually a tree)
 */
typedef struct bptbl_s {
    frame_idx_t  frame;		/**< start or end frame */
    uint8    valid;		/**< For absolute pruning */
    uint8    refcnt;            /**< Reference count (number of successors) */
    int32    wid;		/**< Word index */
    int32    bp;		/**< Back Pointer */
    int32    score;		/**< Score (best among all right contexts) */
    int32    s_idx;		/**< Start of BScoreStack for various right contexts*/
    int32    real_wid;		/**< wid of this or latest predecessor real word */
    int32    prev_real_wid;	/**< wid of second-last real word */
    int16    last_phone;        /**< last phone of this word */
    int16    last2_phone;       /**< next-to-last phone of this word */
} bptbl_t;

/**
 * Segmentation "iterator" for backpointer table results.
 */
typedef struct bptbl_seg_s {
    ps_seg_t base;  /**< Base structure. */
    int32 *bpidx;   /**< Sequence of backpointer IDs. */
    int16 n_bpidx;  /**< Number of backpointer IDs. */
    int16 cur;      /**< Current position in bpidx. */
} bptbl_seg_t;

/*
 * Candidates words for entering their last phones.  Cleared and rebuilt in each
 * frame.
 * NOTE: candidates can only be multi-phone, real dictionary words.
 */
typedef struct lastphn_cand_s {
    int32 wid;
    int32 score;
    int32 bp;
    int32 next;                 /* next candidate starting at the same frame */
} lastphn_cand_t;

/*
 * Since the same instance of a word (i.e., <word,start-frame>) reaches its last
 * phone several times, we can compute its best BP and LM transition score info
 * just the first time and cache it for future occurrences.  Structure for such
 * a cache.
 */
typedef struct {
    int32 sf;                   /* Start frame */
    int32 dscr;                 /* Delta-score upon entering last phone */
    int32 bp;                   /* Best BP */
} last_ltrans_t;

#define CAND_SF_ALLOCSIZE	32
typedef struct {
    int32 bp_ef;
    int32 cand;
} cand_sf_t;

/*
 * Structure for reorganizing the BP table entries in the current frame according
 * to distinct right context ci-phones.  Each entry contains the best BP entry for
 * a given right context.  Each successor word will pick up the correct entry based
 * on its first ci-phone.
 */
typedef struct bestbp_rc_s {
    int32 score;
    int32 path;                 /* BP table index corresponding to this entry */
    int32 lc;                   /* right most ci-phone of above BP entry word */
} bestbp_rc_t;

#define NO_BP		-1

/**
 * Various statistics for profiling.
 */
typedef struct ngram_search_stats_s {
    int32 n_phone_eval;
    int32 n_root_chan_eval;
    int32 n_nonroot_chan_eval;
    int32 n_last_chan_eval;
    int32 n_word_lastchan_eval;
    int32 n_lastphn_cand_utt;
    int32 n_fwdflat_chan;
    int32 n_fwdflat_words;
    int32 n_fwdflat_word_transition;
    int32 n_senone_active_utt;
} ngram_search_stats_t;


/**
 * N-Gram search module structure.
 */
struct ngram_search_s {
    ps_search_t base;
    ngram_model_t *lmset;  /**< Set of language models. */
    hmm_context_t *hmmctx; /**< HMM context. */

    /* Flags to quickly indicate which passes are enabled. */
    uint8 fwdtree;
    uint8 fwdflat;
    uint8 bestpath;

    /* State of procesing. */
    uint8 done;

    /* Allocators */
    listelem_alloc_t *chan_alloc; /**< For chan_t */
    listelem_alloc_t *root_chan_alloc; /**< For root_chan_t */
    listelem_alloc_t *latnode_alloc; /**< For latnode_t */

    /**
     * Search structure of HMM instances.
     *
     * The word triphone sequences (HMM instances) are transformed
     * into tree structures, one tree per unique left triphone in the
     * entire dictionary (actually diphone, since its left context
     * varies dyamically during the search process).  The entire set
     * of trees of channels is allocated once and for all during
     * initialization (since dynamic management of active CHANs is
     * time consuming), with one exception: the last phones of words,
     * that need multiple right context modelling, are not maintained
     * in this static structure since there are too many of them and
     * few are active at any time.  Instead they are maintained as
     * linked lists of CHANs, one list per word, and each CHAN in this
     * set is allocated only on demand and freed if inactive.
     */
    root_chan_t *root_chan;  /**< Roots of search tree. */
    int32 n_root_chan_alloc; /**< Number of root_chan allocated */
    int32 n_root_chan;       /**< Number of valid root_chan */
    int32 n_nonroot_chan;    /**< Number of valid non-root channels */
    int32 max_nonroot_chan;  /**< Maximum possible number of non-root channels */
    root_chan_t *rhmm_1ph;   /**< Root HMMs for single-phone words */

    /**
     * Channels associated with a given word (only used for right
     * contexts, single-phone words in fwdtree search, and word HMMs
     * in fwdflat search).  WARNING: For single-phone words and
     * fwdflat search, this actually contains pointers to root_chan_t,
     * which are allocated using root_chan_alloc.  This is a
     * suboptimal state of affairs.
     */
    chan_t **word_chan;
    bitvec_t *word_active;      /**< array of active flags for all words. */

    /**
     * Each node in the HMM tree structure may point to a set of words
     * whose last phone would follow that node in the tree structure
     * (but is not included in the tree structure for reasons
     * explained above).  The channel node points to one word in this
     * set of words.  The remaining words are linked through
     * homophone_set[].
     * 
     * Single-phone words are not represented in the HMM tree; they
     * are kept in word_chan.
     *
     * Specifically, homophone_set[w] = wid of next word in the same
     * set as w.
     */
    int32 *homophone_set;
    int32 *single_phone_wid; /**< list of single-phone word ids */
    int32 n_1ph_words;       /**< Number single phone words in dict (total) */
    int32 n_1ph_LMwords;     /**< Number single phone dict words also in LM;
                                these come first in single_phone_wid */
    /**
     * Array of active channels for current and next frame.
     *
     * In any frame, only some HMM tree nodes are active.
     * active_chan_list[f mod 2] = list of nonroot channels in the HMM
     * tree active in frame f.
     */
    chan_t ***active_chan_list;
    int32 n_active_chan[2];  /**< Number entries in active_chan_list */
    /**
     * Array of active multi-phone words for current and next frame.
     *
     * Similarly to active_chan_list, active_word_list[f mod 2] = list
     * of word ids for which active channels exist in word_chan in
     * frame f.
     *
     * Statically allocated single-phone words are always active and
     * should not appear in this list.
     */
    int32 **active_word_list;
    int32 n_active_word[2];  /**< Number entries in active_word_list */

    /*
     * FIXME: Document all of these bits.
     */
    lastphn_cand_t *lastphn_cand;
    int32 n_lastphn_cand;
    last_ltrans_t *last_ltrans;      /* one per word */
    int32 cand_sf_alloc;
    cand_sf_t *cand_sf;
    bestbp_rc_t *bestbp_rc;

    bptbl_t *bp_table;       /* Forward pass lattice */
    int32 bpidx;             /* First free BPTable entry */
    int32 bp_table_size;
    int32 *bscore_stack;     /* Score stack for all possible right contexts */
    int32 bss_head;          /* First free BScoreStack entry */
    int32 bscore_stack_size;

    int32 n_frame_alloc; /**< Number of frames allocated in bp_table_idx and friends. */
    int32 n_frame;       /**< Number of frames actually present. */
    int32 *bp_table_idx; /* First BPTable entry for each frame */
    int32 *word_lat_idx; /* BPTable index for any word in current frame;
                            cleared before each frame */

    /*
     * Flat lexicon (2nd pass) search stuff.
     */
    ps_latnode_t **frm_wordlist;   /**< List of active words in each frame. */
    int32 *fwdflat_wordlist;    /**< List of active word IDs for utterance. */
    bitvec_t *expand_word_flag;
    int32 *expand_word_list;
    int32 n_expand_words;
    int32 min_ef_width;
    int32 max_sf_win;
    float32 fwdflat_fwdtree_lw_ratio;

    int32 best_score; /**< Best Viterbi path score. */
    int32 last_phone_best_score; /**< Best Viterbi path score for last phone. */
    int32 renormalized;

    /*
     * DAG (3rd pass) search stuff.
     */
    float32 bestpath_fwdtree_lw_ratio;
    float32 ascale; /**< Acoustic score scale for posterior probabilities. */
    
    ngram_search_stats_t st; /**< Various statistics for profiling. */
    ptmr_t fwdtree_perf;
    ptmr_t fwdflat_perf;
    ptmr_t bestpath_perf;
    int32 n_tot_frame;

    /* A collection of beam widths. */
    int32 beam;
    int32 dynamic_beam;
    int32 pbeam;
    int32 wbeam;
    int32 lpbeam;
    int32 lponlybeam;
    int32 fwdflatbeam;
    int32 fwdflatwbeam;
    int32 fillpen;
    int32 silpen;
    int32 wip;
    int32 nwpen;
    int32 pip;
    int32 maxwpf;
    int32 maxhmmpf;
};
typedef struct ngram_search_s ngram_search_t;

/**
 * Initialize the N-Gram search module.
 */
ps_search_t *ngram_search_init(cmd_ln_t *config,
                               acmod_t *acmod,
                               dict_t *dict,
                               dict2pid_t *d2p);

/**
 * Finalize the N-Gram search module.
 */
void ngram_search_free(ps_search_t *ngs);

/**
 * Record the current frame's index in the backpointer table.
 *
 * @return the current backpointer index.
 */
int ngram_search_mark_bptable(ngram_search_t *ngs, int frame_idx);

/**
 * Enter a word in the backpointer table.
 */
void ngram_search_save_bp(ngram_search_t *ngs, int frame_idx, int32 w,
                          int32 score, int32 path, int32 rc);

/**
 * Allocate last phone channels for all possible right contexts for word w.
 */
void ngram_search_alloc_all_rc(ngram_search_t *ngs, int32 w);

/**
 * Allocate last phone channels for all possible right contexts for word w.
 */
void ngram_search_free_all_rc(ngram_search_t *ngs, int32 w);

/**
 * Find the best word exit for the current frame in the backpointer table.
 *
 * @return the backpointer index of the best word exit.
 */
int ngram_search_find_exit(ngram_search_t *ngs, int frame_idx, int32 *out_best_score, int32 *out_is_final);

/**
 * Backtrace from a given backpointer index to obtain a word hypothesis.
 *
 * @return a <strong>read-only</strong> string with the best hypothesis.
 */
char const *ngram_search_bp_hyp(ngram_search_t *ngs, int bpidx);

/**
 * Compute language and acoustic scores for backpointer table entries.
 */
void ngram_compute_seg_scores(ngram_search_t *ngs, float32 lwf);

/**
 * Construct a word lattice from the current hypothesis.
 */
ps_lattice_t *ngram_search_lattice(ps_search_t *search);

/**
 * Get the exit score for a backpointer entry with a given right context.
 */
int32 ngram_search_exit_score(ngram_search_t *ngs, bptbl_t *pbe, int rcphone);

#endif /* __NGRAM_SEARCH_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file phone_loop_search.h Fast and rough context-independent phoneme loop search.
 */

#include <sphinxbase/err.h>

#include "phone_loop_search.h"

static int phone_loop_search_start(ps_search_t *search);
static int phone_loop_search_step(ps_search_t *search, int frame_idx);
static int phone_loop_search_finish(ps_search_t *search);
static int phone_loop_search_reinit(ps_search_t *search, dict_t *dict, dict2pid_t *d2p);
static void phone_loop_search_free(ps_search_t *search);
static char const *phone_loop_search_hyp(ps_search_t *search, int32 *out_score, int32 *out_is_final);
static int32 phone_loop_search_prob(ps_search_t *search);
static ps_seg_t *phone_loop_search_seg_iter(ps_search_t *search, int32 *out_score);

static ps_searchfuncs_t phone_loop_search_funcs = {
    /* name: */   "phone_loop",
    /* start: */  phone_loop_search_start,
    /* step: */   phone_loop_search_step,
    /* finish: */ phone_loop_search_finish,
    /* reinit: */ phone_loop_search_reinit,
    /* free: */   phone_loop_search_free,
    /* lattice: */  NULL,
    /* hyp: */      phone_loop_search_hyp,
    /* prob: */     phone_loop_search_prob,
    /* seg_iter: */ phone_loop_search_seg_iter,
};

static int
phone_loop_search_reinit(ps_search_t *search, dict_t *dict, dict2pid_t *d2p)
{
    phone_loop_search_t *pls = (phone_loop_search_t *)search;
    cmd_ln_t *config = ps_search_config(search);
    acmod_t *acmod = ps_search_acmod(search);
    int i;

    /* Free old dict2pid, dict, if necessary. */
    ps_search_base_reinit(search, dict, d2p);

    /* Initialize HMM context. */
    if (pls->hmmctx)
        hmm_context_free(pls->hmmctx);
    pls->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
                                   acmod->tmat->tp, NULL, acmod->mdef->sseq);
    if (pls->hmmctx == NULL)
        return -1;

    /* Initialize phone HMMs. */
    if (pls->phones) {
        for (i = 0; i < pls->n_phones; ++i)
            hmm_deinit((hmm_t *)&pls->phones[i]);
        ckd_free(pls->phones);
    }
    pls->n_phones = bin_mdef_n_ciphone(acmod->mdef);
    pls->phones = ckd_calloc(pls->n_phones, sizeof(*pls->phones));
    for (i = 0; i < pls->n_phones; ++i) {
        pls->phones[i].ciphone = i;
        hmm_init(pls->hmmctx, (hmm_t *)&pls->phones[i],
                 FALSE,
                 bin_mdef_pid2ssid(acmod->mdef, i),
                 bin_mdef_pid2tmatid(acmod->mdef, i));
    }
    pls->beam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-pl_beam"));
    pls->pbeam = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-pl_pbeam"));
    pls->pip = logmath_log(acmod->lmath, cmd_ln_float64_r(config, "-pip"));
    E_INFO("State beam %d Phone exit beam %d Insertion penalty %d\n",
           pls->beam, pls->pbeam, pls->pip);

    return 0;
}

ps_search_t *
phone_loop_search_init(cmd_ln_t *config,
		       acmod_t *acmod,
		       dict_t *dict)
{
    phone_loop_search_t *pls;

    /* Allocate and initialize. */
    pls = ckd_calloc(1, sizeof(*pls));
    ps_search_init(ps_search_base(pls), &phone_loop_search_funcs,
                   config, acmod, dict, NULL);
    phone_loop_search_reinit(ps_search_base(pls), ps_search_dict(pls),
                             ps_search_dict2pid(pls));

    return ps_search_base(pls);
}

static void
phone_loop_search_free_renorm(phone_loop_search_t *pls)
{
    gnode_t *gn;
    for (gn = pls->renorm; gn; gn = gnode_next(gn))
        ckd_free(gnode_ptr(gn));
    glist_free(pls->renorm);
    pls->renorm = NULL;
}

static void
phone_loop_search_free(ps_search_t *search)
{
    phone_loop_search_t *pls = (phone_loop_search_t *)search;
    int i;

    ps_search_deinit(search);
    for (i = 0; i < pls->n_phones; ++i)
        hmm_deinit((hmm_t *)&pls->phones[i]);
    phone_loop_search_free_renorm(pls);
    ckd_free(pls->phones);
    hmm_context_free(pls->hmmctx);
    ckd_free(pls);
}

static int
phone_loop_search_start(ps_search_t *search)
{
    phone_loop_search_t *pls = (phone_loop_search_t *)search;
    int i;

    /* Reset and enter all phone HMMs. */
    for (i = 0; i < pls->n_phones; ++i) {
        hmm_t *hmm = (hmm_t *)&pls->phones[i];
        hmm_clear(hmm);
        hmm_enter(hmm, 0, -1, 0);
    }
    phone_loop_search_free_renorm(pls);
    pls->best_score = 0;

    return 0;
}

static void
renormalize_hmms_one(phone_loop_search_t *pls, int frame_idx, int32 norm)
{
    phone_loop_renorm_t *rn = ckd_calloc(1, sizeof(*rn));
    int i;

    pls->renorm = glist_add_ptr(pls->renorm, rn);
    rn->frame_idx = frame_idx;
    rn->norm = norm;

    for (i = 0; i < pls->n_phones; ++i) {
        hmm_normalize((hmm_t *)&pls->phones[i], norm);
    }
}

static int32
evaluate_hmms_one(phone_loop_search_t *pls, int16 const *senscr, int frame_idx)
{
    int32 bs = WORST_SCORE;
    int i;

    hmm_context_set_senscore(pls->hmmctx, senscr);

    for (i = 0; i < pls->n_phones; ++i) {
        hmm_t *hmm = (hmm_t *)&pls->phones[i];
        int32 score;

        if (hmm_frame(hmm) < frame_idx)
            continue;
        score = hmm_vit_eval(hmm);
        if (score BETTER_THAN bs) {
            bs = score;
        }
    }
    pls->best_score = bs;
    return bs;
}

static void
prune_hmms_one(phone_loop_search_t *pls, int frame_idx)
{
    int32 thresh = pls->best_score + pls->beam;
    int nf = frame_idx + 1;
    int i;

    /* Check all phones to see if they remain active in the next frame. */
    for (i = 0; i < pls->n_phones; ++i) {
        hmm_t *hmm = (hmm_t *)&pls->phones[i];

        if (hmm_frame(hmm) < frame_idx)
            continue;
        /* Retain if score better than threshold. */
        if (hmm_bestscore(hmm) BETTER_THAN thresh) {
            hmm_frame(hmm) = nf;
        }
        else
            hmm_clear_scores(hmm);
    }
}

static void
phone_transition_one(phone_loop_search_t *pls, int frame_idx)
{
    int32 thresh = pls->best_score + pls->pbeam;
    int nf = frame_idx + 1;
    int i;

    /* Now transition out of phones whose last states are inside the
     * phone transition beam. */
    for (i = 0; i < pls->n_phones; ++i) {
        hmm_t *hmm = (hmm_t *)&pls->phones[i];
        int32 newphone_score;
        int j;

        if (hmm_frame(hmm) != nf)
            continue;

        newphone_score = hmm_out_score(hmm) + pls->pip;
        if (newphone_score BETTER_THAN thresh) {
            /* Transition into all phones using the usual Viterbi rule. */
            for (j = 0; j < pls->n_phones; ++j) {
                hmm_t *nhmm = (hmm_t *)&pls->phones[j];

                if (hmm_frame(nhmm) < frame_idx
                    || newphone_score BETTER_THAN hmm_in_score(nhmm)) {
                    hmm_enter(nhmm, newphone_score, hmm_out_history(hmm), nf);
                }
            }
        }
    }
}

static int
phone_loop_search_step(ps_search_t *search, int frame_idx)
{
    phone_loop_search_t *pls = (phone_loop_search_t *)search;
    acmod_t *acmod = ps_search_acmod(search);
    int16 const *senscr;
    int i;

    /* All CI senones are active all the time. */
    if (!ps_search_acmod(pls)->compallsen)
        for (i = 0; i < pls->n_phones; ++i)
            acmod_activate_hmm(acmod, (hmm_t *)&pls->phones[i]);

    /* Calculate senone scores for current frame. */
    senscr = acmod_score(acmod, &frame_idx);

    /* Renormalize, if necessary. */
    if (pls->best_score + (2 * pls->beam) WORSE_THAN WORST_SCORE) {
        E_INFO("Renormalizing Scores at frame %d, best score %d\n",
               frame_idx, pls->best_score);
        renormalize_hmms_one(pls, frame_idx, pls->best_score);
    }

    /* Evaluate phone HMMs for current frame. */
    pls->best_score = evaluate_hmms_one(pls, senscr, frame_idx);

    /* Prune phone HMMs. */
    prune_hmms_one(pls, frame_idx);

    /* Do phone transitions. */
    phone_transition_one(pls, frame_idx);

    return 0;
}

static int
phone_loop_search_finish(ps_search_t *search)
{
    /* Actually nothing to do here really. */
    return 0;
}

static char const *
phone_loop_search_hyp(ps_search_t *search, int32 *out_score, int32 *out_is_final)
{
    E_WARN("Hypotheses are not returned from phone loop search");
    return NULL;
}

static int32
phone_loop_search_prob(ps_search_t *search)
{
    /* FIXME: Actually... they ought to be. */
    E_WARN("Posterior probabilities are not returned from phone loop search");
    return 0;
}

static ps_seg_t *
phone_loop_search_seg_iter(ps_search_t *search, int32 *out_score)
{
    E_WARN("Hypotheses are not returned from phone loop search");
    return NULL;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file phone_loop_search.h Fast and rough context-independent
 * phoneme loop search.
 *
 * This exists for the purposes of phoneme lookahead, and thus it
 * actually does not do phoneme recognition (it wouldn't be very
 * accurate anyway).
 */

#ifndef __PHONE_LOOP_SEARCH_H__
#define __PHONE_LOOP_SEARCH_H__

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/ngram_model.h>
#include <sphinxbase/listelem_alloc.h>

/* Local headers. */
#include "pocketsphinx_internal.h"
#include "hmm.h"

/**
 * Phone loop structure.
 */
struct phone_loop_s {
    hmm_t hmm;       /**< Basic HMM structure. */
    int16 ciphone;   /**< Context-independent phone ID. */
    int16 frame;     /**< Last frame this phone was active. */
};
typedef struct phone_loop_s phone_loop_t;

/**
 * Renormalization event.
 */
struct phone_loop_renorm_s {
    int frame_idx;  /**< Frame of renormalization. */
    int32 norm;     /**< Normalization constant. */
};
typedef struct phone_loop_renorm_s phone_loop_renorm_t;

/**
 * Phone loop search structure.
 */
struct phone_loop_search_s {
    ps_search_t base;       /**< Base search structure. */
    hmm_context_t *hmmctx;  /**< HMM context structure. */
    int16 frame;            /**< Current frame being searched. */
    int16 n_phones;         /**< Size of phone array. */
    phone_loop_t *phones;   /**< Array of phone arcs. */

    int32 best_score;       /**< Best Viterbi score in current frame. */
    int32 beam;             /**< HMM pruning beam width. */
    int32 pbeam;            /**< Phone exit pruning beam width. */
    int32 pip;              /**< Phone insertion penalty ("language score"). */
    glist_t renorm;         /**< List of renormalizations. */
};
typedef struct phone_loop_search_s phone_loop_search_t;

ps_search_t *phone_loop_search_init(cmd_ln_t *config,
                                    acmod_t *acmod,
                                    dict_t *dict);

/**
 * Return lookahead heuristic score for a specific phone.
 */
#define phone_loop_search_score(pls,ci) \
    ((pls == NULL) ? 0                                          \
     : (hmm_bestscore(&pls->phones[ci].hmm) - (pls)->best_score))

#endif /* __PHONE_LOOP_SEARCH_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/* System headers. */
#include <stdio.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/err.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/filename.h>
#include <sphinxbase/pio.h>

/* Local headers. */
#include "cmdln_macro.h"
#include "pocketsphinx_internal.h"
#include "ps_lattice_internal.h"
#include "phone_loop_search.h"
#include "fsg_search_internal.h"
#include "ngram_search.h"
#include "ngram_search_fwdtree.h"
#include "ngram_search_fwdflat.h"

static const arg_t ps_args_def[] = {
    POCKETSPHINX_OPTIONS,
    CMDLN_EMPTY_OPTION
};

/* I'm not sure what the portable way to do this is. */
static int
file_exists(const char *path)
{
    FILE *tmp;

    tmp = fopen(path, "rb");
    if (tmp) fclose(tmp);
    return (tmp != NULL);
}

static int
hmmdir_exists(const char *path)
{
    FILE *tmp;
    char *mdef = string_join(path, "/mdef", NULL);

    tmp = fopen(mdef, "rb");
    if (tmp) fclose(tmp);
    ckd_free(mdef);
    return (tmp != NULL);
}

static void
ps_add_file(ps_decoder_t *ps, const char *arg,
            const char *hmmdir, const char *file)
{
    char *tmp = string_join(hmmdir, "/", file, NULL);

    if (cmd_ln_str_r(ps->config, arg) == NULL && file_exists(tmp))
        cmd_ln_set_str_r(ps->config, arg, tmp);
    ckd_free(tmp);
}

static void
ps_init_defaults(ps_decoder_t *ps)
{
    char const *hmmdir, *lmfile, *dictfile;

    /* Disable memory mapping on Blackfin (FIXME: should be uClinux in general). */
#ifdef __ADSPBLACKFIN__
    E_INFO("Will not use mmap() on uClinux/Blackfin.");
    cmd_ln_set_boolean_r(ps->config, "-mmap", FALSE);
#endif

#ifdef MODELDIR
    /* Set default acoustic and language models. */
    hmmdir = cmd_ln_str_r(ps->config, "-hmm");
    lmfile = cmd_ln_str_r(ps->config, "-lm");
    dictfile = cmd_ln_str_r(ps->config, "-dict");
    if (hmmdir == NULL && hmmdir_exists(MODELDIR "/hmm/en_US/hub4wsj_sc_8k")) {
        hmmdir = MODELDIR "/hmm/en_US/hub4wsj_sc_8k";
        cmd_ln_set_str_r(ps->config, "-hmm", hmmdir);
    }
    if (lmfile == NULL && !cmd_ln_str_r(ps->config, "-fsg")
        && !cmd_ln_str_r(ps->config, "-jsgf")
        && file_exists(MODELDIR "/lm/en_US/hub4.5000.DMP")) {
        lmfile = MODELDIR "/lm/en_US/hub4.5000.DMP";
        cmd_ln_set_str_r(ps->config, "-lm", lmfile);
    }
    if (dictfile == NULL && file_exists(MODELDIR "/lm/en_US/cmu07a.dic")) {
        dictfile = MODELDIR "/lm/en_US/cmu07a.dic";
        cmd_ln_set_str_r(ps->config, "-dict", dictfile);
    }

    /* Expand acoustic and language model filenames relative to installation path. */
    if (hmmdir && !path_is_absolute(hmmdir) && !hmmdir_exists(hmmdir)) {
        char *tmphmm = string_join(MODELDIR "/hmm/", hmmdir, NULL);
        if (hmmdir_exists(tmphmm)) {
    	    cmd_ln_set_str_r(ps->config, "-hmm", tmphmm);
    	} else {
    	    E_ERROR("Failed to find mdef file inside the model folder specified with -hmm '%s'\n", hmmdir);
    	}
        ckd_free(tmphmm);
    }
    if (lmfile && !path_is_absolute(lmfile) && !file_exists(lmfile)) {
        char *tmplm = string_join(MODELDIR "/lm/", lmfile, NULL);
        cmd_ln_set_str_r(ps->config, "-lm", tmplm);
        ckd_free(tmplm);
    }
    if (dictfile && !path_is_absolute(dictfile) && !file_exists(dictfile)) {
        char *tmpdict = string_join(MODELDIR "/lm/", dictfile, NULL);
        cmd_ln_set_str_r(ps->config, "-dict", tmpdict);
        ckd_free(tmpdict);
    }
#endif

    /* Get acoustic model filenames and add them to the command-line */
    if ((hmmdir = cmd_ln_str_r(ps->config, "-hmm")) != NULL) {
        ps_add_file(ps, "-mdef", hmmdir, "mdef");
        ps_add_file(ps, "-mean", hmmdir, "means");
        ps_add_file(ps, "-var", hmmdir, "variances");
        ps_add_file(ps, "-tmat", hmmdir, "transition_matrices");
        ps_add_file(ps, "-mixw", hmmdir, "mixture_weights");
        ps_add_file(ps, "-sendump", hmmdir, "sendump");
        ps_add_file(ps, "-fdict", hmmdir, "noisedict");
        ps_add_file(ps, "-lda", hmmdir, "feature_transform");
        ps_add_file(ps, "-featparams", hmmdir, "feat.params");
        ps_add_file(ps, "-senmgau", hmmdir, "senmgau");
    }
}

static void
ps_free_searches(ps_decoder_t *ps)
{
    gnode_t *gn;

    if (ps->searches == NULL)
        return;

    for (gn = ps->searches; gn; gn = gnode_next(gn))
        ps_search_free(gnode_ptr(gn));
    glist_free(ps->searches);
    ps->searches = NULL;
    ps->search = NULL;
}

static ps_search_t *
ps_find_search(ps_decoder_t *ps, char const *name)
{
    gnode_t *gn;

    for (gn = ps->searches; gn; gn = gnode_next(gn)) {
        if (0 == strcmp(ps_search_name(gnode_ptr(gn)), name))
            return (ps_search_t *)gnode_ptr(gn);
    }
    return NULL;
}

int
ps_reinit(ps_decoder_t *ps, cmd_ln_t *config)
{
    char const *lmfile, *lmctl = NULL;

    if (config && config != ps->config) {
        cmd_ln_free_r(ps->config);
        ps->config = cmd_ln_retain(config);
    }

    err_set_debug_level(cmd_ln_int32_r(ps->config, "-debug"));
    ps->mfclogdir = cmd_ln_str_r(ps->config, "-mfclogdir");
    ps->rawlogdir = cmd_ln_str_r(ps->config, "-rawlogdir");
    ps->senlogdir = cmd_ln_str_r(ps->config, "-senlogdir");

    /* Fill in some default arguments. */
    ps_init_defaults(ps);

    /* Free old searches (do this before other reinit) */
    ps_free_searches(ps);

    /* Free old acmod. */
    acmod_free(ps->acmod);
    ps->acmod = NULL;

    /* Free old dictionary (must be done after the two things above) */
    dict_free(ps->dict);
    ps->dict = NULL;
    
    /* Free d2p */
    dict2pid_free(ps->d2p);
    ps->d2p = NULL;

    /* Logmath computation (used in acmod and search) */
    if (ps->lmath == NULL
        || (logmath_get_base(ps->lmath) != 
            (float64)cmd_ln_float32_r(ps->config, "-logbase"))) {
        if (ps->lmath)
            logmath_free(ps->lmath);
        ps->lmath = logmath_init
            ((float64)cmd_ln_float32_r(ps->config, "-logbase"), 0,
             cmd_ln_boolean_r(ps->config, "-bestpath"));
    }

    /* Acoustic model (this is basically everything that
     * uttproc.c, senscr.c, and others used to do) */
    if ((ps->acmod = acmod_init(ps->config, ps->lmath, NULL, NULL)) == NULL)
        return -1;

    if ((ps->pl_window = cmd_ln_int32_r(ps->config, "-pl_window"))) {
        /* Initialize an auxiliary phone loop search, which will run in
         * "parallel" with FSG or N-Gram search. */
        if ((ps->phone_loop = phone_loop_search_init(ps->config,
                                                     ps->acmod, ps->dict)) == NULL)
            return -1;
        ps->searches = glist_add_ptr(ps->searches, ps->phone_loop);
    }

    /* Dictionary and triphone mappings (depends on acmod). */
    /* FIXME: pass config, change arguments, implement LTS, etc. */
    if ((ps->dict = dict_init(ps->config, ps->acmod->mdef)) == NULL)
        return -1;

    /* Determine whether we are starting out in FSG or N-Gram search mode. */
    if (cmd_ln_str_r(ps->config, "-fsg") || cmd_ln_str_r(ps->config, "-jsgf")) {
        ps_search_t *fsgs;

        if ((ps->d2p = dict2pid_build(ps->acmod->mdef, ps->dict)) == NULL)
            return -1;
        if ((fsgs = fsg_search_init(ps->config, ps->acmod, ps->dict, ps->d2p)) == NULL)
            return -1;
        fsgs->pls = ps->phone_loop;
        ps->searches = glist_add_ptr(ps->searches, fsgs);
        ps->search = fsgs;
    }
    else if ((lmfile = cmd_ln_str_r(ps->config, "-lm"))
             || (lmctl = cmd_ln_str_r(ps->config, "-lmctl"))) {
        ps_search_t *ngs;

        /* Make the acmod's feature buffer growable if we are doing two-pass search. */
	if (cmd_ln_boolean_r(ps->config, "-fwdflat")
    	    && cmd_ln_boolean_r(ps->config, "-fwdtree"))
    	    acmod_set_grow(ps->acmod, TRUE);

        if ((ps->d2p = dict2pid_build(ps->acmod->mdef, ps->dict)) == NULL)
            return -1;
        if ((ngs = ngram_search_init(ps->config, ps->acmod, ps->dict, ps->d2p)) == NULL)
            return -1;
        ngs->pls = ps->phone_loop;
        ps->searches = glist_add_ptr(ps->searches, ngs);
        ps->search = ngs;
    }
    /* Otherwise, we will initialize the search whenever the user
     * decides to load an FSG or a language model. */
    else {
        if ((ps->d2p = dict2pid_build(ps->acmod->mdef, ps->dict)) == NULL)
            return -1;
    }

    /* Initialize performance timer. */
    ps->perf.name = "decode";
    ptmr_init(&ps->perf);

    return 0;
}

ps_decoder_t *
ps_init(cmd_ln_t *config)
{
    ps_decoder_t *ps;

    ps = ckd_calloc(1, sizeof(*ps));
    ps->refcount = 1;
    if (ps_reinit(ps, config) < 0) {
        ps_free(ps);
        return NULL;
    }
    return ps;
}

arg_t const *
ps_args(void)
{
    return ps_args_def;
}

ps_decoder_t *
ps_retain(ps_decoder_t *ps)
{
    ++ps->refcount;
    return ps;
}

int
ps_free(ps_decoder_t *ps)
{
    if (ps == NULL)
        return 0;
    if (--ps->refcount > 0)
        return ps->refcount;
    ps_free_searches(ps);
    dict_free(ps->dict);
    dict2pid_free(ps->d2p);
    acmod_free(ps->acmod);
    logmath_free(ps->lmath);
    cmd_ln_free_r(ps->config);
    ckd_free(ps->uttid);
    ckd_free(ps);
    return 0;
}

char const *
ps_get_uttid(ps_decoder_t *ps)
{
    return ps->uttid;
}

cmd_ln_t *
ps_get_config(ps_decoder_t *ps)
{
    return ps->config;
}

logmath_t *
ps_get_logmath(ps_decoder_t *ps)
{
    return ps->lmath;
}

fe_t *
ps_get_fe(ps_decoder_t *ps)
{
    return ps->acmod->fe;
}

feat_t *
ps_get_feat(ps_decoder_t *ps)
{
    return ps->acmod->fcb;
}

ps_mllr_t *
ps_update_mllr(ps_decoder_t *ps, ps_mllr_t *mllr)
{
    return acmod_update_mllr(ps->acmod, mllr);
}

ngram_model_t *
ps_get_lmset(ps_decoder_t *ps)
{
    if (ps->search == NULL
        || 0 != strcmp(ps_search_name(ps->search), "ngram"))
        return NULL;
    return ((ngram_search_t *)ps->search)->lmset;
}

ngram_model_t *
ps_update_lmset(ps_decoder_t *ps, ngram_model_t *lmset)
{
    ngram_search_t *ngs;
    ps_search_t *search;

    /* Look for N-Gram search. */
    search = ps_find_search(ps, "ngram");
    if (search == NULL) {
        /* Initialize N-Gram search. */
        search = ngram_search_init(ps->config, ps->acmod, ps->dict, ps->d2p);
        if (search == NULL)
            return NULL;
        search->pls = ps->phone_loop;
        ps->searches = glist_add_ptr(ps->searches, search);
        ngs = (ngram_search_t *)search;
    }
    else if (lmset != NULL) {
        ngs = (ngram_search_t *)search;
        /* Free any previous lmset if this is a new one. */
        if (ngs->lmset != NULL && ngs->lmset != lmset)
            ngram_model_free(ngs->lmset);
        ngs->lmset = lmset;
        /* Tell N-Gram search to update its view of the world. */
        if (ps_search_reinit(search, ps->dict, ps->d2p) < 0)
            return NULL;
    } else {
	/* Just activate the existing search */
	ngs = (ngram_search_t *)search;
    }
    ps->search = search;
    return ngs->lmset;
}

fsg_set_t *
ps_get_fsgset(ps_decoder_t *ps)
{
    if (ps->search == NULL
        || 0 != strcmp(ps_search_name(ps->search), "fsg"))
        return NULL;
    return (fsg_set_t *)ps->search;
}

fsg_set_t *
ps_update_fsgset(ps_decoder_t *ps)
{
    ps_search_t *search;

    /* Look for FSG search. */
    search = ps_find_search(ps, "fsg");
    if (search == NULL) {
        /* Initialize FSG search. */
        if ((search = fsg_search_init(ps->config,
                                 ps->acmod, ps->dict, ps->d2p)) == NULL) {
            return NULL;
        }
        search->pls = ps->phone_loop;
        ps->searches = glist_add_ptr(ps->searches, search);
    }
    else {
        /* Tell FSG search to update its view of the world. */
        if (ps_search_reinit(search, ps->dict, ps->d2p) < 0)
            return NULL;
    }
    ps->search = search;
    return (fsg_set_t *)search;
}

int
ps_load_dict(ps_decoder_t *ps, char const *dictfile,
             char const *fdictfile, char const *format)
{
    cmd_ln_t *newconfig;
    dict2pid_t *d2p;
    dict_t *dict;
    gnode_t *gn;
    int rv;

    /* Create a new scratch config to load this dict (so existing one
     * won't be affected if it fails) */
    newconfig = cmd_ln_init(NULL, ps_args(), TRUE, NULL);
    cmd_ln_set_boolean_r(newconfig, "-dictcase",
                         cmd_ln_boolean_r(ps->config, "-dictcase"));
    cmd_ln_set_str_r(newconfig, "-dict", dictfile);
    if (fdictfile)
        cmd_ln_set_str_r(newconfig, "-fdict", fdictfile);
    else
        cmd_ln_set_str_r(newconfig, "-fdict",
                         cmd_ln_str_r(ps->config, "-fdict"));

    /* Try to load it. */
    if ((dict = dict_init(newconfig, ps->acmod->mdef)) == NULL) {
        cmd_ln_free_r(newconfig);
        return -1;
    }

    /* Reinit the dict2pid. */
    if ((d2p = dict2pid_build(ps->acmod->mdef, dict)) == NULL) {
        cmd_ln_free_r(newconfig);
        return -1;
    }

    /* Success!  Update the existing config to reflect new dicts and
     * drop everything into place. */
    cmd_ln_free_r(newconfig);
    cmd_ln_set_str_r(ps->config, "-dict", dictfile);
    if (fdictfile)
        cmd_ln_set_str_r(ps->config, "-fdict", fdictfile);
    dict_free(ps->dict);
    ps->dict = dict;
    dict2pid_free(ps->d2p);
    ps->d2p = d2p;

    /* And tell all searches to reconfigure themselves. */
    for (gn = ps->searches; gn; gn = gnode_next(gn)) {
        ps_search_t *search = gnode_ptr(gn);
        if ((rv = ps_search_reinit(search, dict, d2p)) < 0)
            return rv;
    }

    return 0;
}

int
ps_save_dict(ps_decoder_t *ps, char const *dictfile,
             char const *format)
{
    return dict_write(ps->dict, dictfile, format);
}

int
ps_add_word(ps_decoder_t *ps,
            char const *word,
            char const *phones,
            int update)
{
    int32 wid, lmwid;
    ngram_model_t *lmset;
    s3cipid_t *pron;
    char **phonestr, *tmp;
    int np, i, rv;

    /* Parse phones into an array of phone IDs. */
    tmp = ckd_salloc(phones);
    np = str2words(tmp, NULL, 0);
    phonestr = ckd_calloc(np, sizeof(*phonestr));
    str2words(tmp, phonestr, np);
    pron = ckd_calloc(np, sizeof(*pron));
    for (i = 0; i < np; ++i) {
        pron[i] = bin_mdef_ciphone_id(ps->acmod->mdef, phonestr[i]);
        if (pron[i] == -1) {
            E_ERROR("Unknown phone %s in phone string %s\n",
                    phonestr[i], tmp);
            ckd_free(phonestr);
            ckd_free(tmp);
            ckd_free(pron);
            return -1;
        }
    }
    /* No longer needed. */
    ckd_free(phonestr);
    ckd_free(tmp);

    /* Add it to the dictionary. */
    if ((wid = dict_add_word(ps->dict, word, pron, np)) == -1) {
        ckd_free(pron);
        return -1;
    }
    /* No longer needed. */
    ckd_free(pron);

    /* Now we also have to add it to dict2pid. */
    dict2pid_add_word(ps->d2p, wid);

    if ((lmset = ps_get_lmset(ps)) != NULL) {
        /* Add it to the LM set (meaning, the current LM).  In a perfect
         * world, this would result in the same WID, but because of the
         * weird way that word IDs are handled, it doesn't. */
        if ((lmwid = ngram_model_add_word(lmset, word, 1.0))
            == NGRAM_INVALID_WID)
            return -1;
    }
 
    /* Rebuild the widmap and search tree if requested. */
    if (update) {
        if ((rv = ps_search_reinit(ps->search, ps->dict, ps->d2p) < 0))
            return rv;
    }
    return wid;
}

int
ps_decode_raw(ps_decoder_t *ps, FILE *rawfh,
              char const *uttid, long maxsamps)
{
    long total, pos;

    ps_start_utt(ps, uttid);
    /* If this file is seekable or maxsamps is specified, then decode
     * the whole thing at once. */
    if (maxsamps != -1 || (pos = ftell(rawfh)) >= 0) {
        int16 *data;

        if (maxsamps == -1) {
            long endpos;
            fseek(rawfh, 0, SEEK_END);
            endpos = ftell(rawfh);
            fseek(rawfh, pos, SEEK_SET);
            maxsamps = endpos - pos;
        }
        data = ckd_calloc(maxsamps, sizeof(*data));
        total = fread(data, sizeof(*data), maxsamps, rawfh);
        ps_process_raw(ps, data, total, FALSE, TRUE);
        ckd_free(data);
    }
    else {
        /* Otherwise decode it in a stream. */
        total = 0;
        while (!feof(rawfh)) {
            int16 data[256];
            size_t nread;

            nread = fread(data, sizeof(*data), sizeof(data)/sizeof(*data), rawfh);
            ps_process_raw(ps, data, nread, FALSE, FALSE);
            total += nread;
        }
    }
    ps_end_utt(ps);
    return total;
}

int
ps_start_utt(ps_decoder_t *ps, char const *uttid)
{
    int rv;

    if (ps->search == NULL) {
        E_ERROR("No search module is selected, did you forget to "
                "specify a language model or grammar?\n");
        return -1;
    }

    ptmr_reset(&ps->perf);
    ptmr_start(&ps->perf);

    if (uttid) {
        ckd_free(ps->uttid);
        ps->uttid = ckd_salloc(uttid);
    }
    else {
        char nuttid[16];
        ckd_free(ps->uttid);
        sprintf(nuttid, "%09u", ps->uttno);
        ps->uttid = ckd_salloc(nuttid);
        ++ps->uttno;
    }
    /* Remove any residual word lattice and hypothesis. */
    ps_lattice_free(ps->search->dag);
    ps->search->dag = NULL;
    ps->search->last_link = NULL;
    ps->search->post = 0;
    ckd_free(ps->search->hyp_str);
    ps->search->hyp_str = NULL;

    if ((rv = acmod_start_utt(ps->acmod)) < 0)
        return rv;

    /* Start logging features and audio if requested. */
    if (ps->mfclogdir) {
        char *logfn = string_join(ps->mfclogdir, "/",
                                  ps->uttid, ".mfc", NULL);
        FILE *mfcfh;
        E_INFO("Writing MFCC log file: %s\n", logfn);
        if ((mfcfh = fopen(logfn, "wb")) == NULL) {
            E_ERROR_SYSTEM("Failed to open MFCC log file %s", logfn);
            ckd_free(logfn);
            return -1;
        }
        ckd_free(logfn);
        acmod_set_mfcfh(ps->acmod, mfcfh);
    }
    if (ps->rawlogdir) {
        char *logfn = string_join(ps->rawlogdir, "/",
                                  ps->uttid, ".raw", NULL);
        FILE *rawfh;
        E_INFO("Writing raw audio log file: %s\n", logfn);
        if ((rawfh = fopen(logfn, "wb")) == NULL) {
            E_ERROR_SYSTEM("Failed to open raw audio log file %s", logfn);
            ckd_free(logfn);
            return -1;
        }
        ckd_free(logfn);
        acmod_set_rawfh(ps->acmod, rawfh);
    }
    if (ps->senlogdir) {
        char *logfn = string_join(ps->senlogdir, "/",
                                  ps->uttid, ".sen", NULL);
        FILE *senfh;
        E_INFO("Writing senone score log file: %s\n", logfn);
        if ((senfh = fopen(logfn, "wb")) == NULL) {
            E_ERROR_SYSTEM("Failed to open senone score log file %s", logfn);
            ckd_free(logfn);
            return -1;
        }
        ckd_free(logfn);
        acmod_set_senfh(ps->acmod, senfh);
    }

    /* Start auxiliary phone loop search. */
    if (ps->phone_loop)
        ps_search_start(ps->phone_loop);

    return ps_search_start(ps->search);
}

static int
ps_search_forward(ps_decoder_t *ps)
{
    int nfr;

    nfr = 0;
    while (ps->acmod->n_feat_frame > 0) {
        int k;
        if (ps->phone_loop)
            if ((k = ps_search_step(ps->phone_loop, ps->acmod->output_frame)) < 0)
                return k;
        if (ps->acmod->output_frame >= ps->pl_window)
            if ((k = ps_search_step(ps->search,
                                    ps->acmod->output_frame - ps->pl_window)) < 0)
                return k;
        acmod_advance(ps->acmod);
        ++ps->n_frame;
        ++nfr;
    }
    return nfr;
}

int
ps_decode_senscr(ps_decoder_t *ps, FILE *senfh,
                 char const *uttid)
{
    int nfr, n_searchfr;

    ps_start_utt(ps, uttid);
    n_searchfr = 0;
    acmod_set_insenfh(ps->acmod, senfh);
    while ((nfr = acmod_read_scores(ps->acmod)) > 0) {
        if ((nfr = ps_search_forward(ps)) < 0) {
            ps_end_utt(ps);
            return nfr;
        }
        n_searchfr += nfr;
    }
    ps_end_utt(ps);
    acmod_set_insenfh(ps->acmod, NULL);

    return n_searchfr;
}

int
ps_process_raw(ps_decoder_t *ps,
               int16 const *data,
               size_t n_samples,
               int no_search,
               int full_utt)
{
    int n_searchfr = 0;

    if (ps->acmod->state == ACMOD_IDLE) {
	E_ERROR("Failed to process data, utterance is not started. Use start_utt to start it\n");
	return 0;
    }

    if (no_search)
        acmod_set_grow(ps->acmod, TRUE);

    while (n_samples) {
        int nfr;

        /* Process some data into features. */
        if ((nfr = acmod_process_raw(ps->acmod, &data,
                                     &n_samples, full_utt)) < 0)
            return nfr;

        /* Score and search as much data as possible */
        if (no_search)
            continue;
        if ((nfr = ps_search_forward(ps)) < 0)
            return nfr;
        n_searchfr += nfr;
    }

    return n_searchfr;
}

int
ps_process_cep(ps_decoder_t *ps,
               mfcc_t **data,
               int32 n_frames,
               int no_search,
               int full_utt)
{
    int n_searchfr = 0;

    if (no_search)
        acmod_set_grow(ps->acmod, TRUE);

    while (n_frames) {
        int nfr;

        /* Process some data into features. */
        if ((nfr = acmod_process_cep(ps->acmod, &data,
                                     &n_frames, full_utt)) < 0)
            return nfr;

        /* Score and search as much data as possible */
        if (no_search)
            continue;
        if ((nfr = ps_search_forward(ps)) < 0)
            return nfr;
        n_searchfr += nfr;
    }

    return n_searchfr;
}

int
ps_end_utt(ps_decoder_t *ps)
{
    int rv, i;

    acmod_end_utt(ps->acmod);

    /* Search any remaining frames. */
    if ((rv = ps_search_forward(ps)) < 0) {
        ptmr_stop(&ps->perf);
        return rv;
    }
    /* Finish phone loop search. */
    if (ps->phone_loop) {
        if ((rv = ps_search_finish(ps->phone_loop)) < 0) {
            ptmr_stop(&ps->perf);
            return rv;
        }
    }
    /* Search any frames remaining in the lookahead window. */
    for (i = ps->acmod->output_frame - ps->pl_window;
         i < ps->acmod->output_frame; ++i)
        ps_search_step(ps->search, i);
    /* Finish main search. */
    if ((rv = ps_search_finish(ps->search)) < 0) {
        ptmr_stop(&ps->perf);
        return rv;
    }
    ptmr_stop(&ps->perf);

    /* Log a backtrace if requested. */
    if (cmd_ln_boolean_r(ps->config, "-backtrace")) {
        char const *uttid, *hyp;
        ps_seg_t *seg;
        int32 score;

        hyp = ps_get_hyp(ps, &score, &uttid);
        E_INFO("%s: %s (%d)\n", uttid, hyp, score);
        E_INFO_NOFN("%-20s %-5s %-5s %-5s %-10s %-10s %-3s\n",
                    "word", "start", "end", "pprob", "ascr", "lscr", "lback");
        for (seg = ps_seg_iter(ps, &score); seg;
             seg = ps_seg_next(seg)) {
            char const *word;
            int sf, ef;
            int32 post, lscr, ascr, lback;

            word = ps_seg_word(seg);
            ps_seg_frames(seg, &sf, &ef);
            post = ps_seg_prob(seg, &ascr, &lscr, &lback);
            E_INFO_NOFN("%-20s %-5d %-5d %-1.3f %-10d %-10d %-3d\n",
                        word, sf, ef, logmath_exp(ps_get_logmath(ps), post), ascr, lscr, lback);
        }
    }
    return rv;
}

char const *
ps_get_hyp(ps_decoder_t *ps, int32 *out_best_score, char const **out_uttid)
{
    char const *hyp;

    ptmr_start(&ps->perf);
    hyp = ps_search_hyp(ps->search, out_best_score, NULL);
    if (out_uttid)
        *out_uttid = ps->uttid;
    ptmr_stop(&ps->perf);
    return hyp;
}

char const *
ps_get_hyp_final(ps_decoder_t *ps, int32 *out_is_final)
{
    char const *hyp;

    ptmr_start(&ps->perf);
    hyp = ps_search_hyp(ps->search, NULL, out_is_final);
    ptmr_stop(&ps->perf);
    return hyp;
}


int32
ps_get_prob(ps_decoder_t *ps, char const **out_uttid)
{
    int32 prob;

    ptmr_start(&ps->perf);
    prob = ps_search_prob(ps->search);
    if (out_uttid)
        *out_uttid = ps->uttid;
    ptmr_stop(&ps->perf);
    return prob;
}

ps_seg_t *
ps_seg_iter(ps_decoder_t *ps, int32 *out_best_score)
{
    ps_seg_t *itor;

    ptmr_start(&ps->perf);
    itor = ps_search_seg_iter(ps->search, out_best_score);
    ptmr_stop(&ps->perf);
    return itor;
}

ps_seg_t *
ps_seg_next(ps_seg_t *seg)
{
    return ps_search_seg_next(seg);
}

char const *
ps_seg_word(ps_seg_t *seg)
{
    return seg->word;
}

void
ps_seg_frames(ps_seg_t *seg, int *out_sf, int *out_ef)
{
    if (out_sf) *out_sf = seg->sf;
    if (out_ef) *out_ef = seg->ef;
}

int32
ps_seg_prob(ps_seg_t *seg, int32 *out_ascr, int32 *out_lscr, int32 *out_lback)
{
    if (out_ascr) *out_ascr = seg->ascr;
    if (out_lscr) *out_lscr = seg->lscr;
    if (out_lback) *out_lback = seg->lback;
    return seg->prob;
}

void
ps_seg_free(ps_seg_t *seg)
{
    ps_search_seg_free(seg);
}

ps_lattice_t *
ps_get_lattice(ps_decoder_t *ps)
{
    return ps_search_lattice(ps->search);
}

ps_nbest_t *
ps_nbest(ps_decoder_t *ps, int sf, int ef,
         char const *ctx1, char const *ctx2)
{
    ps_lattice_t *dag;
    ngram_model_t *lmset;
    ps_astar_t *nbest;
    float32 lwf;
    int32 w1, w2;

    if (ps->search == NULL)
        return NULL;
    if ((dag = ps_get_lattice(ps)) == NULL)
        return NULL;

    /* FIXME: This is all quite specific to N-Gram search.  Either we
     * should make N-best a method for each search module or it needs
     * to be abstracted to work for N-Gram and FSG. */
    if (0 != strcmp(ps_search_name(ps->search), "ngram")) {
        lmset = NULL;
        lwf = 1.0f;
    }
    else {
        lmset = ((ngram_search_t *)ps->search)->lmset;
        lwf = ((ngram_search_t *)ps->search)->bestpath_fwdtree_lw_ratio;
    }

    w1 = ctx1 ? dict_wordid(ps_search_dict(ps->search), ctx1) : -1;
    w2 = ctx2 ? dict_wordid(ps_search_dict(ps->search), ctx2) : -1;
    nbest = ps_astar_start(dag, lmset, lwf, sf, ef, w1, w2);

    return (ps_nbest_t *)nbest;
}

void
ps_nbest_free(ps_nbest_t *nbest)
{
    ps_astar_finish(nbest);
}

ps_nbest_t *
ps_nbest_next(ps_nbest_t *nbest)
{
    ps_latpath_t *next;

    next = ps_astar_next(nbest);
    if (next == NULL) {
        ps_nbest_free(nbest);
        return NULL;
    }
    return nbest;
}

char const *
ps_nbest_hyp(ps_nbest_t *nbest, int32 *out_score)
{
    assert(nbest != NULL);
    
    if (nbest->top == NULL)
        return NULL;
    if (out_score) *out_score = nbest->top->score;
    return ps_astar_hyp(nbest, nbest->top);
}

ps_seg_t *
ps_nbest_seg(ps_nbest_t *nbest, int32 *out_score)
{
    if (nbest->top == NULL)
        return NULL;
    if (out_score) *out_score = nbest->top->score;
    return ps_astar_seg_iter(nbest, nbest->top, 1.0);
}

int
ps_get_n_frames(ps_decoder_t *ps)
{
    return ps->acmod->output_frame + 1;
}

void
ps_get_utt_time(ps_decoder_t *ps, double *out_nspeech,
                double *out_ncpu, double *out_nwall)
{
    int32 frate;

    frate = cmd_ln_int32_r(ps->config, "-frate");
    *out_nspeech = (double)ps->acmod->output_frame / frate;
    *out_ncpu = ps->perf.t_cpu;
    *out_nwall = ps->perf.t_elapsed;
}

void
ps_get_all_time(ps_decoder_t *ps, double *out_nspeech,
                double *out_ncpu, double *out_nwall)
{
    int32 frate;

    frate = cmd_ln_int32_r(ps->config, "-frate");
    *out_nspeech = (double)ps->n_frame / frate;
    *out_ncpu = ps->perf.t_tot_cpu;
    *out_nwall = ps->perf.t_tot_elapsed;
}

void
ps_search_init(ps_search_t *search, ps_searchfuncs_t *vt,
               cmd_ln_t *config, acmod_t *acmod, dict_t *dict,
               dict2pid_t *d2p)
{
    search->vt = vt;
    search->config = config;
    search->acmod = acmod;
    if (d2p)
        search->d2p = dict2pid_retain(d2p);
    else
        search->d2p = NULL;
    if (dict) {
        search->dict = dict_retain(dict);
        search->start_wid = dict_startwid(dict);
        search->finish_wid = dict_finishwid(dict);
        search->silence_wid = dict_silwid(dict);
        search->n_words = dict_size(dict);
    }
    else {
        search->dict = NULL;
        search->start_wid = search->finish_wid = search->silence_wid = -1;
        search->n_words = 0;
    }
}

void
ps_search_base_reinit(ps_search_t *search, dict_t *dict,
                      dict2pid_t *d2p)
{
    dict_free(search->dict);
    dict2pid_free(search->d2p);
    /* FIXME: _retain() should just return NULL if passed NULL. */
    if (dict) {
        search->dict = dict_retain(dict);
        search->start_wid = dict_startwid(dict);
        search->finish_wid = dict_finishwid(dict);
        search->silence_wid = dict_silwid(dict);
        search->n_words = dict_size(dict);
    }
    else {
        search->dict = NULL;
        search->start_wid = search->finish_wid = search->silence_wid = -1;
        search->n_words = 0;
    }
    if (d2p)
        search->d2p = dict2pid_retain(d2p);
    else
        search->d2p = NULL;
}


void
ps_search_deinit(ps_search_t *search)
{
    /* FIXME: We will have refcounting on acmod, config, etc, at which
     * point we will free them here too. */
    dict_free(search->dict);
    dict2pid_free(search->d2p);
    ckd_free(search->hyp_str);
    ps_lattice_free(search->dag);
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file pocketsphinx_internal.h Internal implementation of
 * PocketSphinx decoder.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __POCKETSPHINX_INTERNAL_H__
#define __POCKETSPHINX_INTERNAL_H__

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/profile.h>

/* Local headers. */
#include "pocketsphinx.h"
#include "acmod.h"
#include "dict.h"
#include "dict2pid.h"

/**
 * Search algorithm structure.
 */
typedef struct ps_search_s ps_search_t;

/**
 * V-table for search algorithm.
 */
typedef struct ps_searchfuncs_s {
    char const *name;

    int (*start)(ps_search_t *search);
    int (*step)(ps_search_t *search, int frame_idx);
    int (*finish)(ps_search_t *search);
    int (*reinit)(ps_search_t *search, dict_t *dict, dict2pid_t *d2p);
    void (*free)(ps_search_t *search);

    ps_lattice_t *(*lattice)(ps_search_t *search);
    char const *(*hyp)(ps_search_t *search, int32 *out_score, int32 *out_is_final);
    int32 (*prob)(ps_search_t *search);
    ps_seg_t *(*seg_iter)(ps_search_t *search, int32 *out_score);
} ps_searchfuncs_t;

/**
 * Base structure for search module.
 */
struct ps_search_s {
    ps_searchfuncs_t *vt;  /**< V-table of search methods. */
    ps_search_t *pls;      /**< Phoneme loop for lookahead. */
    cmd_ln_t *config;      /**< Configuration. */
    acmod_t *acmod;        /**< Acoustic model. */
    dict_t *dict;        /**< Pronunciation dictionary. */
    dict2pid_t *d2p;       /**< Dictionary to senone mappings. */
    char *hyp_str;         /**< Current hypothesis string. */
    ps_lattice_t *dag;	   /**< Current hypothesis word graph. */
    ps_latlink_t *last_link; /**< Final link in best path. */
    int32 post;            /**< Utterance posterior probability. */
    int32 n_words;         /**< Number of words known to search (may
                              be less than in the dictionary) */

    /* Magical word IDs that must exist in the dictionary: */
    int32 start_wid;       /**< Start word ID. */
    int32 silence_wid;     /**< Silence word ID. */
    int32 finish_wid;      /**< Finish word ID. */
};

#define ps_search_base(s) ((ps_search_t *)s)
#define ps_search_config(s) ps_search_base(s)->config
#define ps_search_acmod(s) ps_search_base(s)->acmod
#define ps_search_dict(s) ps_search_base(s)->dict
#define ps_search_dict2pid(s) ps_search_base(s)->d2p
#define ps_search_dag(s) ps_search_base(s)->dag
#define ps_search_last_link(s) ps_search_base(s)->last_link
#define ps_search_post(s) ps_search_base(s)->post
#define ps_search_lookahead(s) ps_search_base(s)->pls
#define ps_search_n_words(s) ps_search_base(s)->n_words

#define ps_search_name(s) ps_search_base(s)->vt->name
#define ps_search_start(s) (*(ps_search_base(s)->vt->start))(s)
#define ps_search_step(s,i) (*(ps_search_base(s)->vt->step))(s,i)
#define ps_search_finish(s) (*(ps_search_base(s)->vt->finish))(s)
#define ps_search_reinit(s,d,d2p) (*(ps_search_base(s)->vt->reinit))(s,d,d2p)
#define ps_search_free(s) (*(ps_search_base(s)->vt->free))(s)
#define ps_search_lattice(s) (*(ps_search_base(s)->vt->lattice))(s)
#define ps_search_hyp(s,sc,final) (*(ps_search_base(s)->vt->hyp))(s,sc,final)
#define ps_search_prob(s) (*(ps_search_base(s)->vt->prob))(s)
#define ps_search_seg_iter(s,sc) (*(ps_search_base(s)->vt->seg_iter))(s,sc)

/* For convenience... */
#define ps_search_silence_wid(s) ps_search_base(s)->silence_wid
#define ps_search_start_wid(s) ps_search_base(s)->start_wid
#define ps_search_finish_wid(s) ps_search_base(s)->finish_wid

/**
 * Initialize base structure.
 */
void ps_search_init(ps_search_t *search, ps_searchfuncs_t *vt,
                    cmd_ln_t *config, acmod_t *acmod, dict_t *dict,
                    dict2pid_t *d2p);

/**
 * Re-initialize base structure with new dictionary.
 */
void ps_search_base_reinit(ps_search_t *search, dict_t *dict,
                           dict2pid_t *d2p);

/**
 * De-initialize base structure.
 */
void ps_search_deinit(ps_search_t *search);

typedef struct ps_segfuncs_s {
    ps_seg_t *(*seg_next)(ps_seg_t *seg);
    void (*seg_free)(ps_seg_t *seg);
} ps_segfuncs_t;

/**
 * Base structure for hypothesis segmentation iterator.
 */
struct ps_seg_s {
    ps_segfuncs_t *vt;     /**< V-table of seg methods */
    ps_search_t *search;   /**< Search object from whence this came */
    char const *word;      /**< Word string (pointer into dictionary hash) */
    frame_idx_t sf;        /**< Start frame. */
    frame_idx_t ef;        /**< End frame. */
    int32 ascr;            /**< Acoustic score. */
    int32 lscr;            /**< Language model score. */
    int32 prob;            /**< Log posterior probability. */
    /* This doesn't need to be 32 bits, so once the scores above are
     * reduced to 16 bits (or less!), this will be too. */
    int32 lback;           /**< Language model backoff. */
    /* Not sure if this should be here at all. */
    float32 lwf;           /**< Language weight factor (for second-pass searches) */
};

#define ps_search_seg_next(seg) (*(seg->vt->seg_next))(seg)
#define ps_search_seg_free(s) (*(seg->vt->seg_free))(seg)


/**
 * Decoder object.
 */
struct ps_decoder_s {
    /* Model parameters and such. */
    cmd_ln_t *config;  /**< Configuration. */
    int refcount;      /**< Reference count. */

    /* Basic units of computation. */
    acmod_t *acmod;    /**< Acoustic model. */
    dict_t *dict;    /**< Pronunciation dictionary. */
    dict2pid_t *d2p;   /**< Dictionary to senone mapping. */
    logmath_t *lmath;  /**< Log math computation. */

    /* Search modules. */
    glist_t searches;        /**< List of search modules. */
    /* TODO: Convert this to a stack of searches each with their own
     * lookahead value. */
    ps_search_t *search;     /**< Currently active search module. */
    ps_search_t *phone_loop; /**< Phone loop search for lookahead. */
    int pl_window;           /**< Window size for phoneme lookahead. */

    /* Utterance-processing related stuff. */
    uint32 uttno;       /**< Utterance counter. */
    char *uttid;        /**< Utterance ID for current utterance. */
    ptmr_t perf;        /**< Performance counter for all of decoding. */
    uint32 n_frame;     /**< Total number of frames processed. */
    char const *mfclogdir; /**< Log directory for MFCC files. */
    char const *rawlogdir; /**< Log directory for audio files. */
    char const *senlogdir; /**< Log directory for senone score files. */
};

#endif /* __POCKETSPHINX_INTERNAL_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ps_alignment.c Multi-level alignment structure
 */

/* System headers. */

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>

/* Local headers. */
#include "ps_alignment.h"

ps_alignment_t *
ps_alignment_init(dict2pid_t *d2p)
{
    ps_alignment_t *al = ckd_calloc(1, sizeof(*al));
    al->d2p = dict2pid_retain(d2p);
    return al;
}

int
ps_alignment_free(ps_alignment_t *al)
{
    if (al == NULL)
        return 0;
    dict2pid_free(al->d2p);
    ckd_free(al->word.seq);
    ckd_free(al->sseq.seq);
    ckd_free(al->state.seq);
    ckd_free(al);
    return 0;
}

#define VECTOR_GROW 10
static void *
vector_grow_one(void *ptr, uint16 *n_alloc, uint16 *n, size_t item_size)
{
    int newsize = *n + 1;
    if (newsize < *n_alloc) {
        *n += 1;
        return ptr;
    }
    newsize += VECTOR_GROW;
    if (newsize > 0xffff)
        return NULL;
    ptr = ckd_realloc(ptr, newsize * item_size);
    *n += 1;
    *n_alloc = newsize;
    return ptr;
}

static ps_alignment_entry_t *
ps_alignment_vector_grow_one(ps_alignment_vector_t *vec)
{
    void *ptr;
    ptr = vector_grow_one(vec->seq, &vec->n_alloc,
                          &vec->n_ent, sizeof(*vec->seq));
    if (ptr == NULL)
        return NULL;
    vec->seq = ptr;
    return vec->seq + vec->n_ent - 1;
}

static void
ps_alignment_vector_empty(ps_alignment_vector_t *vec)
{
    vec->n_ent = 0;
}

int
ps_alignment_add_word(ps_alignment_t *al,
                      int32 wid, int duration)
{
    ps_alignment_entry_t *ent;

    if ((ent = ps_alignment_vector_grow_one(&al->word)) == NULL)
        return 0;
    ent->id.wid = wid;
    if (al->word.n_ent > 1)
        ent->start = ent[-1].start + ent[-1].duration;
    else
        ent->start = 0;
    ent->duration = duration;
    ent->parent = PS_ALIGNMENT_NONE;
    ent->child = PS_ALIGNMENT_NONE;

    return al->word.n_ent;
}

int
ps_alignment_populate(ps_alignment_t *al)
{
    dict2pid_t *d2p;
    dict_t *dict;
    bin_mdef_t *mdef;
    int i, lc;

    /* Clear phone and state sequences. */
    ps_alignment_vector_empty(&al->sseq);
    ps_alignment_vector_empty(&al->state);

    /* For each word, expand to phones/senone sequences. */
    d2p = al->d2p;
    dict = d2p->dict;
    mdef = d2p->mdef;
    lc = bin_mdef_silphone(mdef);
    for (i = 0; i < al->word.n_ent; ++i) {
        ps_alignment_entry_t *went = al->word.seq + i;
        ps_alignment_entry_t *sent;
        int wid = went->id.wid;
        int len = dict_pronlen(dict, wid);
        int j, rc;

        if (i < al->word.n_ent - 1)
            rc = dict_first_phone(dict, al->word.seq[i+1].id.wid);
        else
            rc = bin_mdef_silphone(mdef);

        /* First phone. */
        if ((sent = ps_alignment_vector_grow_one(&al->sseq)) == NULL) {
            E_ERROR("Failed to add phone entry!\n");
            return -1;
        }
        sent->id.pid.cipid = dict_first_phone(dict, wid);
        sent->id.pid.tmatid = bin_mdef_pid2tmatid(mdef, sent->id.pid.cipid);
        sent->start = went->start;
        sent->duration = went->duration;
        sent->parent = i;
        went->child = (uint16)(sent - al->sseq.seq);
        if (len == 1)
            sent->id.pid.ssid
                = dict2pid_lrdiph_rc(d2p, sent->id.pid.cipid, lc, rc);
        else
            sent->id.pid.ssid
                = dict2pid_ldiph_lc(d2p, sent->id.pid.cipid,
                                    dict_second_phone(dict, wid), lc);
        assert(sent->id.pid.ssid != BAD_SSID);

        /* Internal phones. */
        for (j = 1; j < len - 1; ++j) {
            if ((sent = ps_alignment_vector_grow_one(&al->sseq)) == NULL) {
                E_ERROR("Failed to add phone entry!\n");
                return -1;
            }
            sent->id.pid.cipid = dict_pron(dict, wid, j);
            sent->id.pid.tmatid = bin_mdef_pid2tmatid(mdef, sent->id.pid.cipid);
            sent->id.pid.ssid = dict2pid_internal(d2p, wid, j);
            assert(sent->id.pid.ssid != BAD_SSID);
            sent->start = went->start;
            sent->duration = went->duration;
            sent->parent = i;
        }

        /* Last phone. */
        if (j < len) {
            xwdssid_t *rssid;
            assert(j == len - 1);
            if ((sent = ps_alignment_vector_grow_one(&al->sseq)) == NULL) {
                E_ERROR("Failed to add phone entry!\n");
                return -1;
            }
            sent->id.pid.cipid = dict_last_phone(dict, wid);
            sent->id.pid.tmatid = bin_mdef_pid2tmatid(mdef, sent->id.pid.cipid);
            rssid = dict2pid_rssid(d2p, sent->id.pid.cipid,
                                   dict_second_last_phone(dict, wid));
            sent->id.pid.ssid = rssid->ssid[rssid->cimap[rc]];
            assert(sent->id.pid.ssid != BAD_SSID);
            sent->start = went->start;
            sent->duration = went->duration;
            sent->parent = i;
        }
        /* Update lc.  Could just use sent->id.pid.cipid here but that
         * seems needlessly obscure. */
        lc = dict_last_phone(dict, wid);
    }

    /* For each senone sequence, expand to senones.  (we could do this
     * nested above but this makes it more clear and easier to
     * refactor) */
    for (i = 0; i < al->sseq.n_ent; ++i) {
        ps_alignment_entry_t *pent = al->sseq.seq + i;
        ps_alignment_entry_t *sent;
        int j;

        for (j = 0; j < bin_mdef_n_emit_state(mdef); ++j) {
            if ((sent = ps_alignment_vector_grow_one(&al->state)) == NULL) {
                E_ERROR("Failed to add state entry!\n");
                return -1;
            }
            sent->id.senid = bin_mdef_sseq2sen(mdef, pent->id.pid.ssid, j);
            assert(sent->id.senid != BAD_SENID);
            sent->start = pent->start;
            sent->duration = pent->duration;
            sent->parent = i;
            if (j == 0)
                pent->child = (uint16)(sent - al->state.seq);
        }
    }

    return 0;
}

/* FIXME: Somewhat the same as the above function, needs refactoring */
int
ps_alignment_populate_ci(ps_alignment_t *al)
{
    dict2pid_t *d2p;
    dict_t *dict;
    bin_mdef_t *mdef;
    int i;

    /* Clear phone and state sequences. */
    ps_alignment_vector_empty(&al->sseq);
    ps_alignment_vector_empty(&al->state);

    /* For each word, expand to phones/senone sequences. */
    d2p = al->d2p;
    dict = d2p->dict;
    mdef = d2p->mdef;
    for (i = 0; i < al->word.n_ent; ++i) {
        ps_alignment_entry_t *went = al->word.seq + i;
        ps_alignment_entry_t *sent;
        int wid = went->id.wid;
        int len = dict_pronlen(dict, wid);
        int j;

        for (j = 0; j < len; ++j) {
            if ((sent = ps_alignment_vector_grow_one(&al->sseq)) == NULL) {
                E_ERROR("Failed to add phone entry!\n");
                return -1;
            }
            sent->id.pid.cipid = dict_pron(dict, wid, j);
            sent->id.pid.tmatid = bin_mdef_pid2tmatid(mdef, sent->id.pid.cipid);
            sent->id.pid.ssid = bin_mdef_pid2ssid(mdef, sent->id.pid.cipid);
            assert(sent->id.pid.ssid != BAD_SSID);
            sent->start = went->start;
            sent->duration = went->duration;
            sent->parent = i;
        }
    }

    /* For each senone sequence, expand to senones.  (we could do this
     * nested above but this makes it more clear and easier to
     * refactor) */
    for (i = 0; i < al->sseq.n_ent; ++i) {
        ps_alignment_entry_t *pent = al->sseq.seq + i;
        ps_alignment_entry_t *sent;
        int j;

        for (j = 0; j < bin_mdef_n_emit_state(mdef); ++j) {
            if ((sent = ps_alignment_vector_grow_one(&al->state)) == NULL) {
                E_ERROR("Failed to add state entry!\n");
                return -1;
            }
            sent->id.senid = bin_mdef_sseq2sen(mdef, pent->id.pid.ssid, j);
            assert(sent->id.senid != BAD_SENID);
            sent->start = pent->start;
            sent->duration = pent->duration;
            sent->parent = i;
            if (j == 0)
                pent->child = (uint16)(sent - al->state.seq);
        }
    }

    return 0;
}

int
ps_alignment_propagate(ps_alignment_t *al)
{
    ps_alignment_entry_t *last_ent = NULL;
    int i;

    /* Propagate duration up from states to phones. */
    for (i = 0; i < al->state.n_ent; ++i) {
        ps_alignment_entry_t *sent = al->state.seq + i;
        ps_alignment_entry_t *pent = al->sseq.seq + sent->parent;
        if (pent != last_ent) {
            pent->start = sent->start;
            pent->duration = 0;
        }
        pent->duration += sent->duration;
        last_ent = pent;
    }

    /* Propagate duration up from phones to words. */
    last_ent = NULL;
    for (i = 0; i < al->sseq.n_ent; ++i) {
        ps_alignment_entry_t *pent = al->sseq.seq + i;
        ps_alignment_entry_t *went = al->word.seq + pent->parent;
        if (went != last_ent) {
            went->start = pent->start;
            went->duration = 0;
        }
        went->duration += pent->duration;
        last_ent = went;
    }

    return 0;
}

int
ps_alignment_n_words(ps_alignment_t *al)
{
    return (int)al->word.n_ent;
}

int
ps_alignment_n_phones(ps_alignment_t *al)
{
    return (int)al->sseq.n_ent;
}

int
ps_alignment_n_states(ps_alignment_t *al)
{
    return (int)al->state.n_ent;
}

ps_alignment_iter_t *
ps_alignment_words(ps_alignment_t *al)
{
    ps_alignment_iter_t *itor;

    if (al->word.n_ent == 0)
        return NULL;
    itor = ckd_calloc(1, sizeof(*itor));
    itor->al = al;
    itor->vec = &al->word;
    itor->pos = 0;
    return itor;
}

ps_alignment_iter_t *
ps_alignment_phones(ps_alignment_t *al)
{
    ps_alignment_iter_t *itor;

    if (al->sseq.n_ent == 0)
        return NULL;
    itor = ckd_calloc(1, sizeof(*itor));
    itor->al = al;
    itor->vec = &al->sseq;
    itor->pos = 0;
    return itor;
}

ps_alignment_iter_t *
ps_alignment_states(ps_alignment_t *al)
{
    ps_alignment_iter_t *itor;

    if (al->state.n_ent == 0)
        return NULL;
    itor = ckd_calloc(1, sizeof(*itor));
    itor->al = al;
    itor->vec = &al->state;
    itor->pos = 0;
    return itor;
}

ps_alignment_entry_t *
ps_alignment_iter_get(ps_alignment_iter_t *itor)
{
    return itor->vec->seq + itor->pos;
}

int
ps_alignment_iter_free(ps_alignment_iter_t *itor)
{
    ckd_free(itor);
    return 0;
}

ps_alignment_iter_t *
ps_alignment_iter_goto(ps_alignment_iter_t *itor, int pos)
{
    if (itor == NULL)
        return NULL;
    if (pos >= itor->vec->n_ent) {
        ps_alignment_iter_free(itor);
        return NULL;
    }
    itor->pos = pos;
    return itor;
}

ps_alignment_iter_t *
ps_alignment_iter_next(ps_alignment_iter_t *itor)
{
    if (itor == NULL)
        return NULL;
    if (++itor->pos >= itor->vec->n_ent) {
        ps_alignment_iter_free(itor);
        return NULL;
    }
    return itor;
}

ps_alignment_iter_t *
ps_alignment_iter_prev(ps_alignment_iter_t *itor)
{
    if (itor == NULL)
        return NULL;
    if (--itor->pos < 0) {
        ps_alignment_iter_free(itor);
        return NULL;
    }
    return itor;
}

ps_alignment_iter_t *
ps_alignment_iter_up(ps_alignment_iter_t *itor)
{
    ps_alignment_iter_t *itor2;
    if (itor == NULL)
        return NULL;
    if (itor->vec == &itor->al->word)
        return NULL;
    if (itor->vec->seq[itor->pos].parent == PS_ALIGNMENT_NONE)
        return NULL;
    itor2 = ckd_calloc(1, sizeof(*itor2));
    itor2->al = itor->al;
    itor2->pos = itor->vec->seq[itor->pos].parent;
    if (itor->vec == &itor->al->sseq)
        itor2->vec = &itor->al->word;
    else
        itor2->vec = &itor->al->sseq;
    return itor2;
}

ps_alignment_iter_t *
ps_alignment_iter_down(ps_alignment_iter_t *itor)
{
    ps_alignment_iter_t *itor2;
    if (itor == NULL)
        return NULL;
    if (itor->vec == &itor->al->state)
        return NULL;
    if (itor->vec->seq[itor->pos].child == PS_ALIGNMENT_NONE)
        return NULL;
    itor2 = ckd_calloc(1, sizeof(*itor2));
    itor2->al = itor->al;
    itor2->pos = itor->vec->seq[itor->pos].child;
    if (itor->vec == &itor->al->word)
        itor2->vec = &itor->al->sseq;
    else
        itor2->vec = &itor->al->state;
    return itor2;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ps_alignment.h Multi-level alignment structure
 */

#ifndef __PS_ALIGNMENT_H__
#define __PS_ALIGNMENT_H__

/* System headers. */

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>

/* Local headers. */
#include "dict2pid.h"
#include "hmm.h"

#define PS_ALIGNMENT_NONE ((uint16)0xffff)

struct ps_alignment_entry_s {
    union {
        int32 wid;
        struct {
            uint16 ssid;
            uint16 cipid;
            uint16 tmatid;
        } pid;
        uint16 senid;
    } id;
    int16 start;
    int16 duration;
    uint16 parent;
    uint16 child;
};
typedef struct ps_alignment_entry_s ps_alignment_entry_t;

struct ps_alignment_vector_s {
    ps_alignment_entry_t *seq;
    uint16 n_ent, n_alloc;
};
typedef struct ps_alignment_vector_s ps_alignment_vector_t;

struct ps_alignment_s {
    dict2pid_t *d2p;
    ps_alignment_vector_t word;
    ps_alignment_vector_t sseq;
    ps_alignment_vector_t state;
};
typedef struct ps_alignment_s ps_alignment_t;

struct ps_alignment_iter_s {
    ps_alignment_t *al;
    ps_alignment_vector_t *vec;
    int pos;
};
typedef struct ps_alignment_iter_s ps_alignment_iter_t;

/**
 * Create a new, empty alignment.
 */
ps_alignment_t *ps_alignment_init(dict2pid_t *d2p);

/**
 * Release an alignment
 */
int ps_alignment_free(ps_alignment_t *al);

/**
 * Append a word.
 */
int ps_alignment_add_word(ps_alignment_t *al,
                          int32 wid, int duration);

/**
 * Populate lower layers using available word information.
 */
int ps_alignment_populate(ps_alignment_t *al);

/**
 * Populate lower layers using context-independent phones.
 */
int ps_alignment_populate_ci(ps_alignment_t *al);

/**
 * Propagate timing information up from state sequence.
 */
int ps_alignment_propagate(ps_alignment_t *al);

/**
 * Number of words.
 */
int ps_alignment_n_words(ps_alignment_t *al);

/**
 * Number of phones.
 */
int ps_alignment_n_phones(ps_alignment_t *al);

/**
 * Number of states.
 */
int ps_alignment_n_states(ps_alignment_t *al);

/**
 * Iterate over the alignment starting at the first word.
 */
ps_alignment_iter_t *ps_alignment_words(ps_alignment_t *al);

/**
 * Iterate over the alignment starting at the first phone.
 */
ps_alignment_iter_t *ps_alignment_phones(ps_alignment_t *al);

/**
 * Iterate over the alignment starting at the first state.
 */
ps_alignment_iter_t *ps_alignment_states(ps_alignment_t *al);

/**
 * Get the alignment entry pointed to by an iterator.
 */
ps_alignment_entry_t *ps_alignment_iter_get(ps_alignment_iter_t *itor);

/**
 * Move alignment iterator to given index.
 */
ps_alignment_iter_t *ps_alignment_iter_goto(ps_alignment_iter_t *itor, int pos);

/**
 * Move an alignment iterator forward.
 */
ps_alignment_iter_t *ps_alignment_iter_next(ps_alignment_iter_t *itor);

/**
 * Move an alignment iterator back.
 */
ps_alignment_iter_t *ps_alignment_iter_prev(ps_alignment_iter_t *itor);

/**
 * Get a new iterator starting at the parent of the current node.
 */
ps_alignment_iter_t *ps_alignment_iter_up(ps_alignment_iter_t *itor);
/**
 * Get a new iterator starting at the first child of the current node.
 */
ps_alignment_iter_t *ps_alignment_iter_down(ps_alignment_iter_t *itor);

/**
 * Release an iterator before completing all iterations.
 */
int ps_alignment_iter_free(ps_alignment_iter_t *itor);

#endif /* __PS_ALIGNMENT_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ps_lattice.c Word graph search.
 */

/* System headers. */
#include <assert.h>
#include <string.h>
#include <math.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/err.h>
#include <sphinxbase/pio.h>

/* Local headers. */
#include "pocketsphinx_internal.h"
#include "ps_lattice_internal.h"
#include "ngram_search.h"
#include "dict.h"

/*
 * Create a directed link between "from" and "to" nodes, but if a link already exists,
 * choose one with the best ascr.
 */
void
ps_lattice_link(ps_lattice_t *dag, ps_latnode_t *from, ps_latnode_t *to,
                int32 score, int32 ef)
{
    latlink_list_t *fwdlink;

    /* Look for an existing link between "from" and "to" nodes */
    for (fwdlink = from->exits; fwdlink; fwdlink = fwdlink->next)
        if (fwdlink->link->to == to)
            break;

    if (fwdlink == NULL) {
        latlink_list_t *revlink;
        ps_latlink_t *link;

        /* No link between the two nodes; create a new one */
        link = listelem_malloc(dag->latlink_alloc);
        fwdlink = listelem_malloc(dag->latlink_list_alloc);
        revlink = listelem_malloc(dag->latlink_list_alloc);

        link->from = from;
        link->to = to;
        link->ascr = score;
        link->ef = ef;
        link->best_prev = NULL;

        fwdlink->link = revlink->link = link;
        fwdlink->next = from->exits;
        from->exits = fwdlink;
        revlink->next = to->entries;
        to->entries = revlink;
    }
    else {
        /* Link already exists; just retain the best ascr */
        if (score BETTER_THAN fwdlink->link->ascr) {
            fwdlink->link->ascr = score;
            fwdlink->link->ef = ef;
        }
    }           
}

void
ps_lattice_bypass_fillers(ps_lattice_t *dag, int32 silpen, int32 fillpen)
{
    ps_latnode_t *node;
    int32 score;

    /* Bypass filler nodes */
    for (node = dag->nodes; node; node = node->next) {
        latlink_list_t *revlink;
        if (node == dag->end || !dict_filler_word(dag->dict, node->basewid))
            continue;

        /* Replace each link entering filler node with links to all its successors */
        for (revlink = node->entries; revlink; revlink = revlink->next) {
            latlink_list_t *forlink;
            ps_latlink_t *rlink = revlink->link;

            score = (node->basewid == dag->silence) ? silpen : fillpen;
            score += rlink->ascr;
            /*
             * Make links from predecessor of filler (from) to successors of filler.
             * But if successor is a filler, it has already been eliminated since it
             * appears earlier in latnode_list (see build...).  So it can be skipped.
             */
            for (forlink = node->exits; forlink; forlink = forlink->next) {
                ps_latlink_t *flink = forlink->link;
                if (flink->to && rlink->from &&
                    !dict_filler_word(dag->dict, flink->to->basewid)) {
                    ps_lattice_link(dag, rlink->from, flink->to,
                                    score + flink->ascr, flink->ef);
                }
            }
        }
        node->reachable = FALSE;
    }
}

static void
delete_node(ps_lattice_t *dag, ps_latnode_t *node)
{
    latlink_list_t *x, *next_x;

    for (x = node->exits; x; x = next_x) {
        next_x = x->next;
        x->link->from = NULL;
        listelem_free(dag->latlink_list_alloc, x);
    }
    for (x = node->entries; x; x = next_x) {
        next_x = x->next;
        x->link->to = NULL;
        listelem_free(dag->latlink_list_alloc, x);
    }
    listelem_free(dag->latnode_alloc, node);
}


static void
remove_dangling_links(ps_lattice_t *dag, ps_latnode_t *node)
{
    latlink_list_t *x, *prev_x, *next_x;

    prev_x = NULL;
    for (x = node->exits; x; x = next_x) {
        next_x = x->next;
        if (x->link->to == NULL) {
            if (prev_x)
                prev_x->next = next_x;
            else
                node->exits = next_x;
            listelem_free(dag->latlink_alloc, x->link);
            listelem_free(dag->latlink_list_alloc, x);
        }
        else
            prev_x = x;
    }
    prev_x = NULL;
    for (x = node->entries; x; x = next_x) {
        next_x = x->next;
        if (x->link->from == NULL) {
            if (prev_x)
                prev_x->next = next_x;
            else
                node->entries = next_x;
            listelem_free(dag->latlink_alloc, x->link);
            listelem_free(dag->latlink_list_alloc, x);
        }
        else
            prev_x = x;
    }
}

void
ps_lattice_delete_unreachable(ps_lattice_t *dag)
{
    ps_latnode_t *node, *prev_node, *next_node;
    int i;

    /* Remove unreachable nodes from the list of nodes. */
    prev_node = NULL;
    for (node = dag->nodes; node; node = next_node) {
        next_node = node->next;
        if (!node->reachable) {
            if (prev_node)
                prev_node->next = next_node;
            else
                dag->nodes = next_node;
            /* Delete this node and NULLify links to it. */
            delete_node(dag, node);
        }
        else
            prev_node = node;
    }

    /* Remove all links to and from unreachable nodes. */
    i = 0;
    for (node = dag->nodes; node; node = node->next) {
        /* Assign sequence numbers. */
        node->id = i++;

        /* We should obviously not encounter unreachable nodes here! */
        assert(node->reachable);

        /* Remove all links that go nowhere. */
        remove_dangling_links(dag, node);
    }
}

int32
ps_lattice_write(ps_lattice_t *dag, char const *filename)
{
    FILE *fp;
    int32 i;
    ps_latnode_t *d, *initial, *final;

    initial = dag->start;
    final = dag->end;

    E_INFO("Writing lattice file: %s\n", filename);
    if ((fp = fopen(filename, "w")) == NULL) {
        E_ERROR_SYSTEM("Failed to open lattice file '%s' for writing", filename);
        return -1;
    }

    /* Stupid Sphinx-III lattice code expects 'getcwd:' here */
    fprintf(fp, "# getcwd: /this/is/bogus\n");
    fprintf(fp, "# -logbase %e\n", logmath_get_base(dag->lmath));
    fprintf(fp, "#\n");

    fprintf(fp, "Frames %d\n", dag->n_frames);
    fprintf(fp, "#\n");

    for (i = 0, d = dag->nodes; d; d = d->next, i++);
    fprintf(fp,
            "Nodes %d (NODEID WORD STARTFRAME FIRST-ENDFRAME LAST-ENDFRAME)\n",
            i);
    for (i = 0, d = dag->nodes; d; d = d->next, i++) {
        d->id = i;
        fprintf(fp, "%d %s %d %d %d ; %d\n",
                i, dict_wordstr(dag->dict, d->wid),
                d->sf, d->fef, d->lef, d->node_id);
    }
    fprintf(fp, "#\n");

    fprintf(fp, "Initial %d\nFinal %d\n", initial->id, final->id);
    fprintf(fp, "#\n");

    /* Don't bother with this, it's not used by anything. */
    fprintf(fp, "BestSegAscr %d (NODEID ENDFRAME ASCORE)\n",
            0 /* #BPTable entries */ );
    fprintf(fp, "#\n");

    fprintf(fp, "Edges (FROM-NODEID TO-NODEID ASCORE)\n");
    for (d = dag->nodes; d; d = d->next) {
        latlink_list_t *l;
        for (l = d->exits; l; l = l->next) {
            if (l->link->ascr WORSE_THAN WORST_SCORE || l->link->ascr BETTER_THAN 0)
                continue;
            fprintf(fp, "%d %d %d\n",
                    d->id, l->link->to->id, l->link->ascr << SENSCR_SHIFT);
        }
    }
    fprintf(fp, "End\n");
    fclose(fp);

    return 0;
}

int32
ps_lattice_write_htk(ps_lattice_t *dag, char const *filename)
{
    FILE *fp;
    ps_latnode_t *d, *initial, *final;
    int32 j, n_links, n_nodes;

    initial = dag->start;
    final = dag->end;

    E_INFO("Writing lattice file: %s\n", filename);
    if ((fp = fopen(filename, "w")) == NULL) {
        E_ERROR_SYSTEM("Failed to open lattice file '%s' for writing", filename);
        return -1;
    }

    for (n_links = n_nodes = 0, d = dag->nodes; d; d = d->next) {
        latlink_list_t *l;
        if (!d->reachable)
            continue;
        d->id = n_nodes;
        for (l = d->exits; l; l = l->next) {
            if (l->link->to == NULL || !l->link->to->reachable)
                continue;
            if (l->link->ascr WORSE_THAN WORST_SCORE || l->link->ascr BETTER_THAN 0)
                continue;
            ++n_links;
        }
        ++n_nodes;
    }

    fprintf(fp, "# Lattice generated by PocketSphinx\n");
    fprintf(fp, "#\n# Header\n#\n");
    fprintf(fp, "VERSION=1.0\n");
    fprintf(fp, "start=%d\n", initial->id);
    fprintf(fp, "end=%d\n", final->id);
    fprintf(fp, "#\n");

    fprintf(fp, "N=%d\tL=%d\n", n_nodes, n_links);
    fprintf(fp, "#\n# Node definitions\n#\n");
    for (d = dag->nodes; d; d = d->next) {
        char const *word = dict_wordstr(dag->dict, d->wid);
        char const *c = strrchr(word, '(');
        int altpron = 1;
        if (!d->reachable)
            continue;
        if (c)
            altpron = atoi(c + 1);
        word = dict_basestr(dag->dict, d->wid);
        if (d->wid == dict_startwid(dag->dict))
            word = "!SENT_START";
        else if (d->wid == dict_finishwid(dag->dict))
            word = "!SENT_END";
        else if (dict_filler_word(dag->dict, d->wid))
            word = "!NULL";
        fprintf(fp, "I=%d\tt=%.2f\tW=%s\tv=%d\n",
                d->id, (double)d->sf / dag->frate,
                word, altpron);
    }
    fprintf(fp, "#\n# Link definitions\n#\n");
    for (j = 0, d = dag->nodes; d; d = d->next) {
        latlink_list_t *l;
        if (!d->reachable)
            continue;
        for (l = d->exits; l; l = l->next) {
            if (l->link->to == NULL || !l->link->to->reachable)
                continue;
            if (l->link->ascr WORSE_THAN WORST_SCORE || l->link->ascr BETTER_THAN 0)
                continue;
            fprintf(fp, "J=%d\tS=%d\tE=%d\ta=%f\tp=%g\n", j++,
                    d->id, l->link->to->id,
                    logmath_log_to_ln(dag->lmath, l->link->ascr << SENSCR_SHIFT),
                    logmath_exp(dag->lmath, l->link->alpha + l->link->beta - dag->norm));
        }
    }
    fclose(fp);

    return 0;
}

/* Read parameter from a lattice file*/
static int
dag_param_read(lineiter_t *li, char *param)
{
    int32 n;

    while ((li = lineiter_next(li)) != NULL) {
        char *c;

        /* Ignore comments. */
        if (li->buf[0] == '#')
            continue;

        /* Find the first space. */
        c = strchr(li->buf, ' ');
        if (c == NULL) continue;

        /* Check that the first field equals param and that there's a number after it. */
        if (strncmp(li->buf, param, strlen(param)) == 0
            && sscanf(c + 1, "%d", &n) == 1)
            return n;
    }
    return -1;
}

/* Mark every node that has a path to the argument dagnode as "reachable". */
static void
dag_mark_reachable(ps_latnode_t * d)
{
    latlink_list_t *l;

    d->reachable = 1;
    for (l = d->entries; l; l = l->next)
        if (l->link->from && !l->link->from->reachable)
            dag_mark_reachable(l->link->from);
}

ps_lattice_t *
ps_lattice_read(ps_decoder_t *ps,
                char const *file)
{
    FILE *fp;
    int32 ispipe;
    lineiter_t *line;
    float64 lb;
    float32 logratio;
    ps_latnode_t *tail;
    ps_latnode_t **darray;
    ps_lattice_t *dag;
    int i, k, n_nodes;
    int32 pip, silpen, fillpen;

    dag = ckd_calloc(1, sizeof(*dag));

    if (ps) {
        dag->search = ps->search;
        dag->dict = dict_retain(ps->dict);
        dag->lmath = logmath_retain(ps->lmath);
        dag->frate = cmd_ln_int32_r(dag->search->config, "-frate");
    }
    else {
        dag->dict = dict_init(NULL, NULL);
        dag->lmath = logmath_init(1.0001, 0, FALSE);
        dag->frate = 100;
    }
    dag->silence = dict_silwid(dag->dict);
    dag->latnode_alloc = listelem_alloc_init(sizeof(ps_latnode_t));
    dag->latlink_alloc = listelem_alloc_init(sizeof(ps_latlink_t));
    dag->latlink_list_alloc = listelem_alloc_init(sizeof(latlink_list_t));
    dag->refcount = 1;

    tail = NULL;
    darray = NULL;

    E_INFO("Reading DAG file: %s\n", file);
    if ((fp = fopen_compchk(file, &ispipe)) == NULL) {
        E_ERROR_SYSTEM("Failed to open DAG file '%s' for reading", file);
        return NULL;
    }
    line = lineiter_start(fp);

    /* Read and verify logbase (ONE BIG HACK!!) */
    if (line == NULL) {
        E_ERROR("Premature EOF(%s)\n", file);
        goto load_error;
    }
    if (strncmp(line->buf, "# getcwd: ", 10) != 0) {
        E_ERROR("%s does not begin with '# getcwd: '\n%s", file, line->buf);
        goto load_error;
    }
    if ((line = lineiter_next(line)) == NULL) {
        E_ERROR("Premature EOF(%s)\n", file);
        goto load_error;
    }
    if ((strncmp(line->buf, "# -logbase ", 11) != 0)
        || (sscanf(line->buf + 11, "%lf", &lb) != 1)) {
        E_WARN("%s: Cannot find -logbase in header\n", file);
        lb = 1.0001;
    }
    logratio = 1.0f;
    if (dag->lmath == NULL)
        dag->lmath = logmath_init(lb, 0, TRUE);
    else {
        float32 pb = logmath_get_base(dag->lmath);
        if (fabs(lb - pb) >= 0.0001) {
            E_WARN("Inconsistent logbases: %f vs %f: will compensate\n", lb, pb);
            logratio = (float32)(log(lb) / log(pb));
            E_INFO("Lattice log ratio: %f\n", logratio);
        }
    }
    /* Read Frames parameter */
    dag->n_frames = dag_param_read(line, "Frames");
    if (dag->n_frames <= 0) {
        E_ERROR("Frames parameter missing or invalid\n");
        goto load_error;
    }
    /* Read Nodes parameter */
    n_nodes = dag_param_read(line, "Nodes");
    if (n_nodes <= 0) {
        E_ERROR("Nodes parameter missing or invalid\n");
        goto load_error;
    }

    /* Read nodes */
    darray = ckd_calloc(n_nodes, sizeof(*darray));
    for (i = 0; i < n_nodes; i++) {
        ps_latnode_t *d;
        int32 w;
        int seqid, sf, fef, lef;
        char wd[256];

        if ((line = lineiter_next(line)) == NULL) {
            E_ERROR("Premature EOF while loading Nodes(%s)\n", file);
            goto load_error;
        }

        if ((k =
             sscanf(line->buf, "%d %255s %d %d %d", &seqid, wd, &sf, &fef,
                    &lef)) != 5) {
            E_ERROR("Cannot parse line: %s, value of count %d\n", line->buf, k);
            goto load_error;
        }

        w = dict_wordid(dag->dict, wd);
        if (w < 0) {
            if (dag->search == NULL) {
                char *ww = ckd_salloc(wd);
                if (dict_word2basestr(ww) != -1) {
                    if (dict_wordid(dag->dict, ww) == BAD_S3WID)
                        dict_add_word(dag->dict, ww, NULL, 0);
                }
                ckd_free(ww);
                w = dict_add_word(dag->dict, wd, NULL, 0);
            }
            if (w < 0) {
                E_ERROR("Unknown word in line: %s\n", line->buf);
                goto load_error;
            }
        }

        if (seqid != i) {
            E_ERROR("Seqno error: %s\n", line->buf);
            goto load_error;
        }

        d = listelem_malloc(dag->latnode_alloc);
        darray[i] = d;
        d->wid = w;
        d->basewid = dict_basewid(dag->dict, w);
        d->id = seqid;
        d->sf = sf;
        d->fef = fef;
        d->lef = lef;
        d->reachable = 0;
        d->exits = d->entries = NULL;
        d->next = NULL;

        if (!dag->nodes)
            dag->nodes = d;
        else
            tail->next = d;
        tail = d;
    }

    /* Read initial node ID */
    k = dag_param_read(line, "Initial");
    if ((k < 0) || (k >= n_nodes)) {
        E_ERROR("Initial node parameter missing or invalid\n");
        goto load_error;
    }
    dag->start = darray[k];

    /* Read final node ID */
    k = dag_param_read(line, "Final");
    if ((k < 0) || (k >= n_nodes)) {
        E_ERROR("Final node parameter missing or invalid\n");
        goto load_error;
    }
    dag->end = darray[k];

    /* Read bestsegscore entries and ignore them. */
    if ((k = dag_param_read(line, "BestSegAscr")) < 0) {
        E_ERROR("BestSegAscr parameter missing\n");
        goto load_error;
    }
    for (i = 0; i < k; i++) {
        if ((line = lineiter_next(line)) == NULL) {
            E_ERROR("Premature EOF while (%s) ignoring BestSegAscr\n",
                    line);
            goto load_error;
        }
    }

    /* Read in edges. */
    while ((line = lineiter_next(line)) != NULL) {
        if (line->buf[0] == '#')
            continue;
        if (0 == strncmp(line->buf, "Edges", 5))
            break;
    }
    if (line == NULL) {
        E_ERROR("Edges missing\n");
        goto load_error;
    }
    while ((line = lineiter_next(line)) != NULL) {
        int from, to, ascr;
        ps_latnode_t *pd, *d;

        if (sscanf(line->buf, "%d %d %d", &from, &to, &ascr) != 3)
            break;
        if (ascr WORSE_THAN WORST_SCORE)
            continue;
        pd = darray[from];
        d = darray[to];
        if (logratio != 1.0f)
            ascr = (int32)(ascr * logratio);
        ps_lattice_link(dag, pd, d, ascr, d->sf - 1);
    }
    if (strcmp(line->buf, "End\n") != 0) {
        E_ERROR("Terminating 'End' missing\n");
        goto load_error;
    }
    lineiter_free(line);
    fclose_comp(fp, ispipe);
    ckd_free(darray);

    /* Minor hack: If the final node is a filler word and not </s>,
     * then set its base word ID to </s>, so that the language model
     * scores won't be screwed up. */
    if (dict_filler_word(dag->dict, dag->end->wid))
        dag->end->basewid = dag->search
            ? ps_search_finish_wid(dag->search)
            : dict_wordid(dag->dict, S3_FINISH_WORD);

    /* Mark reachable from dag->end */
    dag_mark_reachable(dag->end);

    /* Free nodes unreachable from dag->end and their links */
    ps_lattice_delete_unreachable(dag);

    if (ps) {
        /* Build links around silence and filler words, since they do
         * not exist in the language model.  FIXME: This is
         * potentially buggy, as we already do this before outputing
         * lattices. */
        pip = logmath_log(dag->lmath, cmd_ln_float32_r(ps->config, "-pip"));
        silpen = pip + logmath_log(dag->lmath,
                                   cmd_ln_float32_r(ps->config, "-silprob"));
        fillpen = pip + logmath_log(dag->lmath,
                                    cmd_ln_float32_r(ps->config, "-fillprob"));
        ps_lattice_bypass_fillers(dag, silpen, fillpen);
    }

    return dag;

  load_error:
    E_ERROR("Failed to load %s\n", file);
    lineiter_free(line);
    if (fp) fclose_comp(fp, ispipe);
    ckd_free(darray);
    return NULL;
}

int
ps_lattice_n_frames(ps_lattice_t *dag)
{
    return dag->n_frames;
}

ps_lattice_t *
ps_lattice_init_search(ps_search_t *search, int n_frame)
{
    ps_lattice_t *dag;

    dag = ckd_calloc(1, sizeof(*dag));
    dag->search = search;
    dag->dict = dict_retain(search->dict);
    dag->lmath = logmath_retain(search->acmod->lmath);
    dag->frate = cmd_ln_int32_r(dag->search->config, "-frate");
    dag->silence = dict_silwid(dag->dict);
    dag->n_frames = n_frame;
    dag->latnode_alloc = listelem_alloc_init(sizeof(ps_latnode_t));
    dag->latlink_alloc = listelem_alloc_init(sizeof(ps_latlink_t));
    dag->latlink_list_alloc = listelem_alloc_init(sizeof(latlink_list_t));
    dag->refcount = 1;
    return dag;
}

ps_lattice_t *
ps_lattice_retain(ps_lattice_t *dag)
{
    ++dag->refcount;
    return dag;
}

int
ps_lattice_free(ps_lattice_t *dag)
{
    if (dag == NULL)
        return 0;
    if (--dag->refcount > 0)
        return dag->refcount;
    logmath_free(dag->lmath);
    dict_free(dag->dict);
    listelem_alloc_free(dag->latnode_alloc);
    listelem_alloc_free(dag->latlink_alloc);
    listelem_alloc_free(dag->latlink_list_alloc);    
    ckd_free(dag->hyp_str);
    ckd_free(dag);
    return 0;
}

logmath_t *
ps_lattice_get_logmath(ps_lattice_t *dag)
{
    return dag->lmath;
}

ps_latnode_iter_t *
ps_latnode_iter(ps_lattice_t *dag)
{
    return dag->nodes;
}

ps_latnode_iter_t *
ps_latnode_iter_next(ps_latnode_iter_t *itor)
{
    return itor->next;
}

void
ps_latnode_iter_free(ps_latnode_iter_t *itor)
{
    /* Do absolutely nothing. */
}

ps_latnode_t *
ps_latnode_iter_node(ps_latnode_iter_t *itor)
{
    return itor;
}

int
ps_latnode_times(ps_latnode_t *node, int16 *out_fef, int16 *out_lef)
{
    if (out_fef) *out_fef = (int16)node->fef;
    if (out_lef) *out_lef = (int16)node->lef;
    return node->sf;
}

char const *
ps_latnode_word(ps_lattice_t *dag, ps_latnode_t *node)
{
    return dict_wordstr(dag->dict, node->wid);
}

char const *
ps_latnode_baseword(ps_lattice_t *dag, ps_latnode_t *node)
{
    return dict_wordstr(dag->dict, node->basewid);
}

int32
ps_latnode_prob(ps_lattice_t *dag, ps_latnode_t *node,
                ps_latlink_t **out_link)
{
    latlink_list_t *links;
    int32 bestpost = logmath_get_zero(dag->lmath);

    for (links = node->exits; links; links = links->next) {
        int32 post = links->link->alpha + links->link->beta - dag->norm;
        if (post > bestpost) {
            if (out_link) *out_link = links->link;
            bestpost = post;
        }
    }
    return bestpost;
}

ps_latlink_iter_t *
ps_latnode_exits(ps_latnode_t *node)
{
    return node->exits;
}

ps_latlink_iter_t *
ps_latnode_entries(ps_latnode_t *node)
{
    return node->entries;
}

ps_latlink_iter_t *
ps_latlink_iter_next(ps_latlink_iter_t *itor)
{
    return itor->next;
}

void
ps_latlink_iter_free(ps_latlink_iter_t *itor)
{
    /* Do absolutely nothing. */
}

ps_latlink_t *
ps_latlink_iter_link(ps_latlink_iter_t *itor)
{
    return itor->link;
}

int
ps_latlink_times(ps_latlink_t *link, int16 *out_sf)
{
    if (out_sf) {
        if (link->from) {
            *out_sf = link->from->sf;
        }
        else {
            *out_sf = 0;
        }
    }
    return link->ef;
}

ps_latnode_t *
ps_latlink_nodes(ps_latlink_t *link, ps_latnode_t **out_src)
{
    if (out_src) *out_src = link->from;
    return link->to;
}

char const *
ps_latlink_word(ps_lattice_t *dag, ps_latlink_t *link)
{
    if (link->from == NULL)
        return NULL;
    return dict_wordstr(dag->dict, link->from->wid);
}

char const *
ps_latlink_baseword(ps_lattice_t *dag, ps_latlink_t *link)
{
    if (link->from == NULL)
        return NULL;
    return dict_wordstr(dag->dict, link->from->basewid);
}

ps_latlink_t *
ps_latlink_pred(ps_latlink_t *link)
{
    return link->best_prev;
}

int32
ps_latlink_prob(ps_lattice_t *dag, ps_latlink_t *link, int32 *out_ascr)
{
    int32 post = link->alpha + link->beta - dag->norm;
    if (out_ascr) *out_ascr = link->ascr << SENSCR_SHIFT;
    return post;
}

char const *
ps_lattice_hyp(ps_lattice_t *dag, ps_latlink_t *link)
{
    ps_latlink_t *l;
    size_t len;
    char *c;

    /* Backtrace once to get hypothesis length. */
    len = 0;
    /* FIXME: There may not be a search, but actually there should be a dict. */
    if (dict_real_word(dag->dict, link->to->basewid)) {
	char *wstr = dict_wordstr(dag->dict, link->to->basewid);
	if (wstr != NULL)
	    len += strlen(wstr) + 1;
    }
    for (l = link; l; l = l->best_prev) {
        if (dict_real_word(dag->dict, l->from->basewid)) {
    	    char *wstr = dict_wordstr(dag->dict, l->from->basewid);
            if (wstr != NULL)
        	len += strlen(wstr) + 1;
        }
    }

    /* Backtrace again to construct hypothesis string. */
    ckd_free(dag->hyp_str);
    dag->hyp_str = ckd_calloc(1, len+1); /* extra one incase the hyp is empty */
    c = dag->hyp_str + len - 1;
    if (dict_real_word(dag->dict, link->to->basewid)) {
	char *wstr = dict_wordstr(dag->dict, link->to->basewid);
	if (wstr != NULL) {
    	    len = strlen(wstr);
	    c -= len;
    	    memcpy(c, wstr, len);
    	    if (c > dag->hyp_str) {
        	--c;
        	*c = ' ';
	    }
        }
    }
    for (l = link; l; l = l->best_prev) {
        if (dict_real_word(dag->dict, l->from->basewid)) {
    	    char *wstr = dict_wordstr(dag->dict, l->from->basewid);
    	    if (wstr != NULL) {
	        len = strlen(wstr);            
    		c -= len;
    		memcpy(c, wstr, len);
        	if (c > dag->hyp_str) {
            	    --c;
            	    *c = ' ';
        	}
    	    }
        }
    }

    return dag->hyp_str;
}

static void
ps_lattice_compute_lscr(ps_seg_t *seg, ps_latlink_t *link, int to)
{
    ngram_model_t *lmset;

    /* Language model score is included in the link score for FSG
     * search.  FIXME: Of course, this is sort of a hack :( */
    if (0 != strcmp(ps_search_name(seg->search), "ngram")) {
        seg->lback = 1; /* Unigram... */
        seg->lscr = 0;
        return;
    }
        
    lmset = ((ngram_search_t *)seg->search)->lmset;

    if (link->best_prev == NULL) {
        if (to) /* Sentence has only two words. */
            seg->lscr = ngram_bg_score(lmset, link->to->basewid,
                                       link->from->basewid, &seg->lback)
                >> SENSCR_SHIFT;
        else {/* This is the start symbol, its lscr is always 0. */
            seg->lscr = 0;
            seg->lback = 1;
        }
    }
    else {
        /* Find the two predecessor words. */
        if (to) {
            seg->lscr = ngram_tg_score(lmset, link->to->basewid,
                                       link->from->basewid,
                                       link->best_prev->from->basewid,
                                       &seg->lback) >> SENSCR_SHIFT;
        }
        else {
            if (link->best_prev->best_prev)
                seg->lscr = ngram_tg_score(lmset, link->from->basewid,
                                           link->best_prev->from->basewid,
                                           link->best_prev->best_prev->from->basewid,
                                           &seg->lback) >> SENSCR_SHIFT;
            else
                seg->lscr = ngram_bg_score(lmset, link->from->basewid,
                                           link->best_prev->from->basewid,
                                           &seg->lback) >> SENSCR_SHIFT;
        }
    }
}

static void
ps_lattice_link2itor(ps_seg_t *seg, ps_latlink_t *link, int to)
{
    dag_seg_t *itor = (dag_seg_t *)seg;
    ps_latnode_t *node;

    if (to) {
        node = link->to;
        seg->ef = node->lef;
        seg->prob = 0; /* norm + beta - norm */
    }
    else {
        latlink_list_t *x;
        ps_latnode_t *n;
        logmath_t *lmath = ps_search_acmod(seg->search)->lmath;

        node = link->from;
        seg->ef = link->ef;
        seg->prob = link->alpha + link->beta - itor->norm;
        /* Sum over all exits for this word and any alternate
           pronunciations at the same frame. */
        for (n = node; n; n = n->alt) {
            for (x = n->exits; x; x = x->next) {
                if (x->link == link)
                    continue;
                seg->prob = logmath_add(lmath, seg->prob,
                                        x->link->alpha + x->link->beta - itor->norm);
            }
        }
    }
    seg->word = dict_wordstr(ps_search_dict(seg->search), node->wid);
    seg->sf = node->sf;
    seg->ascr = link->ascr << SENSCR_SHIFT;
    /* Compute language model score from best predecessors. */
    ps_lattice_compute_lscr(seg, link, to);
}

static void
ps_lattice_seg_free(ps_seg_t *seg)
{
    dag_seg_t *itor = (dag_seg_t *)seg;
    
    ckd_free(itor->links);
    ckd_free(itor);
}

static ps_seg_t *
ps_lattice_seg_next(ps_seg_t *seg)
{
    dag_seg_t *itor = (dag_seg_t *)seg;

    ++itor->cur;
    if (itor->cur == itor->n_links + 1) {
        ps_lattice_seg_free(seg);
        return NULL;
    }
    else if (itor->cur == itor->n_links) {
        /* Re-use the last link but with the "to" node. */
        ps_lattice_link2itor(seg, itor->links[itor->cur - 1], TRUE);
    }
    else {
        ps_lattice_link2itor(seg, itor->links[itor->cur], FALSE);
    }

    return seg;
}

static ps_segfuncs_t ps_lattice_segfuncs = {
    /* seg_next */ ps_lattice_seg_next,
    /* seg_free */ ps_lattice_seg_free
};

ps_seg_t *
ps_lattice_seg_iter(ps_lattice_t *dag, ps_latlink_t *link, float32 lwf)
{
    dag_seg_t *itor;
    ps_latlink_t *l;
    int cur;

    /* Calling this an "iterator" is a bit of a misnomer since we have
     * to get the entire backtrace in order to produce it.
     */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &ps_lattice_segfuncs;
    itor->base.search = dag->search;
    itor->base.lwf = lwf;
    itor->n_links = 0;
    itor->norm = dag->norm;

    for (l = link; l; l = l->best_prev) {
        ++itor->n_links;
    }
    if (itor->n_links == 0) {
        ckd_free(itor);
        return NULL;
    }

    itor->links = ckd_calloc(itor->n_links, sizeof(*itor->links));
    cur = itor->n_links - 1;
    for (l = link; l; l = l->best_prev) {
        itor->links[cur] = l;
        --cur;
    }

    ps_lattice_link2itor((ps_seg_t *)itor, itor->links[0], FALSE);
    return (ps_seg_t *)itor;
}

latlink_list_t *
latlink_list_new(ps_lattice_t *dag, ps_latlink_t *link, latlink_list_t *next)
{
    latlink_list_t *ll;

    ll = listelem_malloc(dag->latlink_list_alloc);
    ll->link = link;
    ll->next = next;

    return ll;
}

void
ps_lattice_pushq(ps_lattice_t *dag, ps_latlink_t *link)
{
    if (dag->q_head == NULL)
        dag->q_head = dag->q_tail = latlink_list_new(dag, link, NULL);
    else {
        dag->q_tail->next = latlink_list_new(dag, link, NULL);
        dag->q_tail = dag->q_tail->next;
    }

}

ps_latlink_t *
ps_lattice_popq(ps_lattice_t *dag)
{
    latlink_list_t *x;
    ps_latlink_t *link;

    if (dag->q_head == NULL)
        return NULL;
    link = dag->q_head->link;
    x = dag->q_head->next;
    listelem_free(dag->latlink_list_alloc, dag->q_head);
    dag->q_head = x;
    if (dag->q_head == NULL)
        dag->q_tail = NULL;
    return link;
}

void
ps_lattice_delq(ps_lattice_t *dag)
{
    while (ps_lattice_popq(dag)) {
        /* Do nothing. */
    }
}

ps_latlink_t *
ps_lattice_traverse_edges(ps_lattice_t *dag, ps_latnode_t *start, ps_latnode_t *end)
{
    ps_latnode_t *node;
    latlink_list_t *x;

    /* Cancel any unfinished traversal. */
    ps_lattice_delq(dag);

    /* Initialize node fanin counts and path scores. */
    for (node = dag->nodes; node; node = node->next)
        node->info.fanin = 0;
    for (node = dag->nodes; node; node = node->next) {
        for (x = node->exits; x; x = x->next)
            (x->link->to->info.fanin)++;
    }

    /* Initialize agenda with all exits from start. */
    if (start == NULL) start = dag->start;
    for (x = start->exits; x; x = x->next)
        ps_lattice_pushq(dag, x->link);

    /* Pull the first edge off the queue. */
    return ps_lattice_traverse_next(dag, end);
}

ps_latlink_t *
ps_lattice_traverse_next(ps_lattice_t *dag, ps_latnode_t *end)
{
    ps_latlink_t *next;

    next = ps_lattice_popq(dag);
    if (next == NULL)
        return NULL;

    /* Decrease fanin count for destination node and expand outgoing
     * edges if all incoming edges have been seen. */
    --next->to->info.fanin;
    if (next->to->info.fanin == 0) {
        latlink_list_t *x;

        if (end == NULL) end = dag->end;
        if (next->to == end) {
            /* If we have traversed all links entering the end node,
             * clear the queue, causing future calls to this function
             * to return NULL. */
            ps_lattice_delq(dag);
            return next;
        }

        /* Extend all outgoing edges. */
        for (x = next->to->exits; x; x = x->next)
            ps_lattice_pushq(dag, x->link);
    }
    return next;
}

ps_latlink_t *
ps_lattice_reverse_edges(ps_lattice_t *dag, ps_latnode_t *start, ps_latnode_t *end)
{
    ps_latnode_t *node;
    latlink_list_t *x;

    /* Cancel any unfinished traversal. */
    ps_lattice_delq(dag);

    /* Initialize node fanout counts and path scores. */
    for (node = dag->nodes; node; node = node->next) {
        node->info.fanin = 0;
        for (x = node->exits; x; x = x->next)
            ++node->info.fanin;
    }

    /* Initialize agenda with all entries from end. */
    if (end == NULL) end = dag->end;
    for (x = end->entries; x; x = x->next)
        ps_lattice_pushq(dag, x->link);

    /* Pull the first edge off the queue. */
    return ps_lattice_reverse_next(dag, start);
}

ps_latlink_t *
ps_lattice_reverse_next(ps_lattice_t *dag, ps_latnode_t *start)
{
    ps_latlink_t *next;

    next = ps_lattice_popq(dag);
    if (next == NULL)
        return NULL;

    /* Decrease fanout count for source node and expand incoming
     * edges if all incoming edges have been seen. */
    --next->from->info.fanin;
    if (next->from->info.fanin == 0) {
        latlink_list_t *x;

        if (start == NULL) start = dag->start;
        if (next->from == start) {
            /* If we have traversed all links entering the start node,
             * clear the queue, causing future calls to this function
             * to return NULL. */
            ps_lattice_delq(dag);
            return next;
        }

        /* Extend all outgoing edges. */
        for (x = next->from->entries; x; x = x->next)
            ps_lattice_pushq(dag, x->link);
    }
    return next;
}

/*
 * Find the best score from dag->start to end point of any link and
 * use it to update links further down the path.  This is like
 * single-source shortest path search, except that it is done over
 * edges rather than nodes, which allows us to do exact trigram scoring.
 *
 * Helpfully enough, we get half of the posterior probability
 * calculation for free that way too.  (interesting research topic: is
 * there a reliable Viterbi analogue to word-level Forward-Backward
 * like there is for state-level?  Or, is it just lattice density?)
 */
ps_latlink_t *
ps_lattice_bestpath(ps_lattice_t *dag, ngram_model_t *lmset,
                    float32 lwf, float32 ascale)
{
    ps_search_t *search;
    ps_latnode_t *node;
    ps_latlink_t *link;
    ps_latlink_t *bestend;
    latlink_list_t *x;
    logmath_t *lmath;
    int32 bestescr;

    search = dag->search;
    lmath = dag->lmath;

    /* Initialize path scores for all links exiting dag->start, and
     * set all other scores to the minimum.  Also initialize alphas to
     * log-zero. */
    for (node = dag->nodes; node; node = node->next) {
        for (x = node->exits; x; x = x->next) {
            x->link->path_scr = MAX_NEG_INT32;
            x->link->alpha = logmath_get_zero(lmath);
        }
    }
    for (x = dag->start->exits; x; x = x->next) {
        int32 n_used;

        /* Ignore filler words. */
        if (dict_filler_word(ps_search_dict(search), x->link->to->basewid)
            && x->link->to != dag->end)
            continue;

        /* Best path points to dag->start, obviously. */
        if (lmset)
            x->link->path_scr = x->link->ascr +
                (ngram_bg_score(lmset, x->link->to->basewid,
                                ps_search_start_wid(search), &n_used) 
                 >> SENSCR_SHIFT)
                 * lwf;
        else
            x->link->path_scr = x->link->ascr;
        x->link->best_prev = NULL;
        /* No predecessors for start links. */
        x->link->alpha = 0;
    }

    /* Traverse the edges in the graph, updating path scores. */
    for (link = ps_lattice_traverse_edges(dag, NULL, NULL);
         link; link = ps_lattice_traverse_next(dag, NULL)) {
        int32 bprob, n_used;

        /* Skip filler nodes in traversal. */
        if (dict_filler_word(ps_search_dict(search), link->from->basewid) && link->from != dag->start)
            continue;
        if (dict_filler_word(ps_search_dict(search), link->to->basewid) && link->to != dag->end)
            continue;

        /* Sanity check, we should not be traversing edges that
         * weren't previously updated, otherwise nasty overflows will result. */
        assert(link->path_scr != MAX_NEG_INT32);

        /* Calculate common bigram probability for all alphas. */
        if (lmset)
            bprob = ngram_ng_prob(lmset,
                                  link->to->basewid,
                                  &link->from->basewid, 1, &n_used);
        else
            bprob = 0;
        /* Add in this link's acoustic score, which was a constant
           factor in previous computations (if any). */
        link->alpha += (link->ascr << SENSCR_SHIFT) * ascale;

        /* Update scores for all paths exiting link->to. */
        for (x = link->to->exits; x; x = x->next) {
            int32 tscore, score;

            /* Skip links to filler words in update. */
            if (dict_filler_word(ps_search_dict(search), x->link->to->basewid)
                && x->link->to != dag->end)
                continue;

            /* Update alpha with sum of previous alphas. */
            x->link->alpha = logmath_add(lmath, x->link->alpha, link->alpha + bprob);
            /* Calculate trigram score for bestpath. */
            if (lmset)
                tscore = (ngram_tg_score(lmset, x->link->to->basewid,
                                        link->to->basewid,
                                        link->from->basewid, &n_used) >> SENSCR_SHIFT)
                    * lwf;
            else
                tscore = 0;
            /* Update link score with maximum link score. */
            score = link->path_scr + tscore + x->link->ascr;
            if (score BETTER_THAN x->link->path_scr) {
                x->link->path_scr = score;
                x->link->best_prev = link;
            }
        }
    }

    /* Find best link entering final node, and calculate normalizer
     * for posterior probabilities. */
    bestend = NULL;
    bestescr = MAX_NEG_INT32;

    /* Normalizer is the alpha for the imaginary link exiting the
       final node. */
    dag->norm = logmath_get_zero(lmath);
    for (x = dag->end->entries; x; x = x->next) {
        int32 bprob, n_used;

        if (dict_filler_word(ps_search_dict(search), x->link->from->basewid))
            continue;
        if (lmset)
            bprob = ngram_ng_prob(lmset,
                                  x->link->to->basewid,
                                  &x->link->from->basewid, 1, &n_used);
        else
            bprob = 0;
        dag->norm = logmath_add(lmath, dag->norm, x->link->alpha + bprob);
        if (x->link->path_scr BETTER_THAN bestescr) {
            bestescr = x->link->path_scr;
            bestend = x->link;
        }
    }
    /* FIXME: floating point... */
    dag->norm += (int32)(dag->final_node_ascr << SENSCR_SHIFT) * ascale;

    E_INFO("Normalizer P(O) = alpha(%s:%d:%d) = %d\n",
           dict_wordstr(dag->search->dict, dag->end->wid),
           dag->end->sf, dag->end->lef,
           dag->norm);
    return bestend;
}

static int32
ps_lattice_joint(ps_lattice_t *dag, ps_latlink_t *link, float32 ascale)
{
    ngram_model_t *lmset;
    int32 jprob;

    /* Sort of a hack... */
    if (dag->search && 0 == strcmp(ps_search_name(dag->search), "ngram"))
        lmset = ((ngram_search_t *)dag->search)->lmset;
    else
        lmset = NULL;

    jprob = (dag->final_node_ascr << SENSCR_SHIFT) * ascale;
    while (link) {
        if (lmset) {
            int lback;
            /* Compute unscaled language model probability.  Note that
               this is actually not the language model probability
               that corresponds to this link, but that is okay,
               because we are just taking the sum over all links in
               the best path. */
            jprob += ngram_ng_prob(lmset, link->to->basewid,
                                   &link->from->basewid, 1, &lback);
        }
        /* If there is no language model, we assume that the language
           model probability (such as it is) has been included in the
           link score. */
        jprob += (link->ascr << SENSCR_SHIFT) * ascale;
        link = link->best_prev;
    }

    E_INFO("Joint P(O,S) = %d P(S|O) = %d\n", jprob, jprob - dag->norm);
    return jprob;
}

int32
ps_lattice_posterior(ps_lattice_t *dag, ngram_model_t *lmset,
                     float32 ascale)
{
    ps_search_t *search;
    logmath_t *lmath;
    ps_latnode_t *node;
    ps_latlink_t *link;
    latlink_list_t *x;
    ps_latlink_t *bestend;
    int32 bestescr;

    search = dag->search;
    lmath = dag->lmath;

    /* Reset all betas to zero. */
    for (node = dag->nodes; node; node = node->next) {
        for (x = node->exits; x; x = x->next) {
            x->link->beta = logmath_get_zero(lmath);
        }
    }

    bestend = NULL;
    bestescr = MAX_NEG_INT32;
    /* Accumulate backward probabilities for all links. */
    for (link = ps_lattice_reverse_edges(dag, NULL, NULL);
         link; link = ps_lattice_reverse_next(dag, NULL)) {
        int32 bprob, n_used;

        /* Skip filler nodes in traversal. */
        if (dict_filler_word(ps_search_dict(search), link->from->basewid) && link->from != dag->start)
            continue;
        if (dict_filler_word(ps_search_dict(search), link->to->basewid) && link->to != dag->end)
            continue;

        /* Calculate LM probability. */
        if (lmset)
            bprob = ngram_ng_prob(lmset, link->to->basewid,
                                  &link->from->basewid, 1, &n_used);
        else
            bprob = 0;

        if (link->to == dag->end) {
            /* Track the best path - we will backtrace in order to
               calculate the unscaled joint probability for sentence
               posterior. */
            if (link->path_scr BETTER_THAN bestescr) {
                bestescr = link->path_scr;
                bestend = link;
            }
            /* Imaginary exit link from final node has beta = 1.0 */
            link->beta = bprob + (dag->final_node_ascr << SENSCR_SHIFT) * ascale;
        }
        else {
            /* Update beta from all outgoing betas. */
            for (x = link->to->exits; x; x = x->next) {
                if (dict_filler_word(ps_search_dict(search), x->link->to->basewid) && x->link->to != dag->end)
                    continue;
                link->beta = logmath_add(lmath, link->beta,
                                         x->link->beta + bprob
                                         + (x->link->ascr << SENSCR_SHIFT) * ascale);
            }
        }
    }

    /* Return P(S|O) = P(O,S)/P(O) */
    return ps_lattice_joint(dag, bestend, ascale) - dag->norm;
}

int32
ps_lattice_posterior_prune(ps_lattice_t *dag, int32 beam)
{
    ps_latlink_t *link;
    int npruned = 0;

    for (link = ps_lattice_traverse_edges(dag, dag->start, dag->end);
         link; link = ps_lattice_traverse_next(dag, dag->end)) {
        link->from->reachable = FALSE;
        if (link->alpha + link->beta - dag->norm < beam) {
            latlink_list_t *x, *tmp, *next;
            tmp = NULL;
            for (x = link->from->exits; x; x = next) {
                next = x->next;
                if (x->link == link) {
                    listelem_free(dag->latlink_list_alloc, x);
                }
                else {
                    x->next = tmp;
                    tmp = x;
                }
            }
            link->from->exits = tmp;
            tmp = NULL;
            for (x = link->to->entries; x; x = next) {
                next = x->next;
                if (x->link == link) {
                    listelem_free(dag->latlink_list_alloc, x);
                }
                else {
                    x->next = tmp;
                    tmp = x;
                }
            }
            link->to->entries = tmp;
            listelem_free(dag->latlink_alloc, link);
            ++npruned;
        }
    }
    dag_mark_reachable(dag->end);
    ps_lattice_delete_unreachable(dag);
    return npruned;
}


/* Parameters to prune n-best alternatives search */
#define MAX_PATHS	500     /* Max allowed active paths at any time */
#define MAX_HYP_TRIES	10000

/*
 * For each node in any path between from and end of utt, find the
 * best score from "from".sf to end of utt.  (NOTE: Uses bigram probs;
 * this is an estimate of the best score from "from".)  (NOTE #2: yes,
 * this is the "heuristic score" used in A* search)
 */
static int32
best_rem_score(ps_astar_t *nbest, ps_latnode_t * from)
{
    latlink_list_t *x;
    int32 bestscore, score;

    if (from->info.rem_score <= 0)
        return (from->info.rem_score);

    /* Best score from "from" to end of utt not known; compute from successors */
    bestscore = WORST_SCORE;
    for (x = from->exits; x; x = x->next) {
        int32 n_used;

        score = best_rem_score(nbest, x->link->to);
        score += x->link->ascr;
        if (nbest->lmset)
            score += (ngram_bg_score(nbest->lmset, x->link->to->basewid,
                                     from->basewid, &n_used) >> SENSCR_SHIFT)
                      * nbest->lwf;
        if (score BETTER_THAN bestscore)
            bestscore = score;
    }
    from->info.rem_score = bestscore;

    return bestscore;
}

/*
 * Insert newpath in sorted (by path score) list of paths.  But if newpath is
 * too far down the list, drop it (FIXME: necessary?)
 * total_score = path score (newpath) + rem_score to end of utt.
 */
static void
path_insert(ps_astar_t *nbest, ps_latpath_t *newpath, int32 total_score)
{
    ps_latpath_t *prev, *p;
    int32 i;

    prev = NULL;
    for (i = 0, p = nbest->path_list; (i < MAX_PATHS) && p; p = p->next, i++) {
        if ((p->score + p->node->info.rem_score) < total_score)
            break;
        prev = p;
    }

    /* newpath should be inserted between prev and p */
    if (i < MAX_PATHS) {
        /* Insert new partial hyp */
        newpath->next = p;
        if (!prev)
            nbest->path_list = newpath;
        else
            prev->next = newpath;
        if (!p)
            nbest->path_tail = newpath;

        nbest->n_path++;
        nbest->n_hyp_insert++;
        nbest->insert_depth += i;
    }
    else {
        /* newpath score too low; reject it and also prune paths beyond MAX_PATHS */
        nbest->path_tail = prev;
        prev->next = NULL;
        nbest->n_path = MAX_PATHS;
        listelem_free(nbest->latpath_alloc, newpath);

        nbest->n_hyp_reject++;
        for (; p; p = newpath) {
            newpath = p->next;
            listelem_free(nbest->latpath_alloc, p);
            nbest->n_hyp_reject++;
        }
    }
}

/* Find all possible extensions to given partial path */
static void
path_extend(ps_astar_t *nbest, ps_latpath_t * path)
{
    latlink_list_t *x;
    ps_latpath_t *newpath;
    int32 total_score, tail_score;

    /* Consider all successors of path->node */
    for (x = path->node->exits; x; x = x->next) {
        int32 n_used;

        /* Skip successor if no path from it reaches the final node */
        if (x->link->to->info.rem_score <= WORST_SCORE)
            continue;

        /* Create path extension and compute exact score for this extension */
        newpath = listelem_malloc(nbest->latpath_alloc);
        newpath->node = x->link->to;
        newpath->parent = path;
        newpath->score = path->score + x->link->ascr;
        if (nbest->lmset) {
            if (path->parent) {
                newpath->score += nbest->lwf
                    * (ngram_tg_score(nbest->lmset, newpath->node->basewid,
                                      path->node->basewid,
                                      path->parent->node->basewid, &n_used)
                       >> SENSCR_SHIFT);
            }
            else 
                newpath->score += nbest->lwf
                    * (ngram_bg_score(nbest->lmset, newpath->node->basewid,
                                      path->node->basewid, &n_used)
                       >> SENSCR_SHIFT);
        }

        /* Insert new partial path hypothesis into sorted path_list */
        nbest->n_hyp_tried++;
        total_score = newpath->score + newpath->node->info.rem_score;

        /* First see if hyp would be worse than the worst */
        if (nbest->n_path >= MAX_PATHS) {
            tail_score =
                nbest->path_tail->score
                + nbest->path_tail->node->info.rem_score;
            if (total_score < tail_score) {
                listelem_free(nbest->latpath_alloc, newpath);
                nbest->n_hyp_reject++;
                continue;
            }
        }

        path_insert(nbest, newpath, total_score);
    }
}

ps_astar_t *
ps_astar_start(ps_lattice_t *dag,
                  ngram_model_t *lmset,
                  float32 lwf,
                  int sf, int ef,
                  int w1, int w2)
{
    ps_astar_t *nbest;
    ps_latnode_t *node;

    nbest = ckd_calloc(1, sizeof(*nbest));
    nbest->dag = dag;
    nbest->lmset = lmset;
    nbest->lwf = lwf;
    nbest->sf = sf;
    if (ef < 0)
        nbest->ef = dag->n_frames + 1;
    else
        nbest->ef = ef;
    nbest->w1 = w1;
    nbest->w2 = w2;
    nbest->latpath_alloc = listelem_alloc_init(sizeof(ps_latpath_t));

    /* Initialize rem_score (A* heuristic) to default values */
    for (node = dag->nodes; node; node = node->next) {
        if (node == dag->end)
            node->info.rem_score = 0;
        else if (node->exits == NULL)
            node->info.rem_score = WORST_SCORE;
        else
            node->info.rem_score = 1;   /* +ve => unknown value */
    }

    /* Create initial partial hypotheses list consisting of nodes starting at sf */
    nbest->path_list = nbest->path_tail = NULL;
    for (node = dag->nodes; node; node = node->next) {
        if (node->sf == sf) {
            ps_latpath_t *path;
            int32 n_used;

            best_rem_score(nbest, node);
            path = listelem_malloc(nbest->latpath_alloc);
            path->node = node;
            path->parent = NULL;
            if (nbest->lmset)
                path->score = nbest->lwf *
                    ((w1 < 0)
                    ? ngram_bg_score(nbest->lmset, node->basewid, w2, &n_used)
                    : ngram_tg_score(nbest->lmset, node->basewid, w2, w1, &n_used));
            else
                path->score = 0;
            path->score >>= SENSCR_SHIFT;
            path_insert(nbest, path, path->score + node->info.rem_score);
        }
    }

    return nbest;
}

ps_latpath_t *
ps_astar_next(ps_astar_t *nbest)
{
    ps_lattice_t *dag;

    dag = nbest->dag;

    /* Pop the top (best) partial hypothesis */
    while ((nbest->top = nbest->path_list) != NULL) {
        nbest->path_list = nbest->path_list->next;
        if (nbest->top == nbest->path_tail)
            nbest->path_tail = NULL;
        nbest->n_path--;

        /* Complete hypothesis? */
        if ((nbest->top->node->sf >= nbest->ef)
            || ((nbest->top->node == dag->end) &&
                (nbest->ef > dag->end->sf))) {
            /* FIXME: Verify that it is non-empty.  Also we may want
             * to verify that it is actually distinct from other
             * paths, since often this is not the case*/
            return nbest->top;
        }
        else {
            if (nbest->top->node->fef < nbest->ef)
                path_extend(nbest, nbest->top);
        }
    }

    /* Did not find any more paths to extend. */
    return NULL;
}

char const *
ps_astar_hyp(ps_astar_t *nbest, ps_latpath_t *path)
{
    ps_search_t *search;
    ps_latpath_t *p;
    size_t len;
    char *c;
    char *hyp;

    search = nbest->dag->search;

    /* Backtrace once to get hypothesis length. */
    len = 0;
    for (p = path; p; p = p->parent) {
        if (dict_real_word(ps_search_dict(search), p->node->basewid)) {
    	    char *wstr = dict_wordstr(ps_search_dict(search), p->node->basewid);
    	    if (wstr != NULL)
    	        len += strlen(wstr) + 1;
        }
    }

    if (len == 0) {
	return NULL;
    }

    /* Backtrace again to construct hypothesis string. */
    hyp = ckd_calloc(1, len);
    c = hyp + len - 1;
    for (p = path; p; p = p->parent) {
        if (dict_real_word(ps_search_dict(search), p->node->basewid)) {
    	    char *wstr = dict_wordstr(ps_search_dict(search), p->node->basewid);
    	    if (wstr != NULL) {
	        len = strlen(wstr);
    		c -= len;
        	memcpy(c, wstr, len);
    		if (c > hyp) {
            	    --c;
        	    *c = ' ';
    		}
    	    }
        }
    }

    nbest->hyps = glist_add_ptr(nbest->hyps, hyp);
    return hyp;
}

static void
ps_astar_node2itor(astar_seg_t *itor)
{
    ps_seg_t *seg = (ps_seg_t *)itor;
    ps_latnode_t *node;

    assert(itor->cur < itor->n_nodes);
    node = itor->nodes[itor->cur];
    if (itor->cur == itor->n_nodes - 1)
        seg->ef = node->lef;
    else
        seg->ef = itor->nodes[itor->cur + 1]->sf - 1;
    seg->word = dict_wordstr(ps_search_dict(seg->search), node->wid);
    seg->sf = node->sf;
    seg->prob = 0; /* FIXME: implement forward-backward */
}

static void
ps_astar_seg_free(ps_seg_t *seg)
{
    astar_seg_t *itor = (astar_seg_t *)seg;
    ckd_free(itor->nodes);
    ckd_free(itor);
}

static ps_seg_t *
ps_astar_seg_next(ps_seg_t *seg)
{
    astar_seg_t *itor = (astar_seg_t *)seg;

    ++itor->cur;
    if (itor->cur == itor->n_nodes) {
        ps_astar_seg_free(seg);
        return NULL;
    }
    else {
        ps_astar_node2itor(itor);
    }

    return seg;
}

static ps_segfuncs_t ps_astar_segfuncs = {
    /* seg_next */ ps_astar_seg_next,
    /* seg_free */ ps_astar_seg_free
};

ps_seg_t *
ps_astar_seg_iter(ps_astar_t *astar, ps_latpath_t *path, float32 lwf)
{
    astar_seg_t *itor;
    ps_latpath_t *p;
    int cur;

    /* Backtrace and make an iterator, this should look familiar by now. */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &ps_astar_segfuncs;
    itor->base.search = astar->dag->search;
    itor->base.lwf = lwf;
    itor->n_nodes = itor->cur = 0;
    for (p = path; p; p = p->parent) {
        ++itor->n_nodes;
    }
    itor->nodes = ckd_calloc(itor->n_nodes, sizeof(*itor->nodes));
    cur = itor->n_nodes - 1;
    for (p = path; p; p = p->parent) {
        itor->nodes[cur] = p->node;
        --cur;
    }

    ps_astar_node2itor(itor);
    return (ps_seg_t *)itor;
}

void
ps_astar_finish(ps_astar_t *nbest)
{
    gnode_t *gn;

    /* Free all hyps. */
    for (gn = nbest->hyps; gn; gn = gnode_next(gn)) {
        ckd_free(gnode_ptr(gn));
    }
    glist_free(nbest->hyps);
    /* Free all paths. */
    listelem_alloc_free(nbest->latpath_alloc);
    /* Free the Henge. */
    ckd_free(nbest);
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ps_lattice_internal.h Word graph search implementation
 */

#ifndef __PS_LATTICE_INTERNAL_H__
#define __PS_LATTICE_INTERNAL_H__

/**
 * Linked list of DAG link pointers.
 *
 * Because the same link structure is used for forward and reverse
 * links, as well as for the agenda used in bestpath search, we can't
 * store the list pointer inside latlink_t.  We could use glist_t
 * here, but it wastes 4 bytes per entry on 32-bit machines.
 */
typedef struct latlink_list_s {
    ps_latlink_t *link;
    struct latlink_list_s *next;
} latlink_list_t;

/**
 * Word graph structure used in bestpath/nbest search.
 */
struct ps_lattice_s {
    int refcount;      /**< Reference count. */

    logmath_t *lmath;    /**< Log-math object. */
    ps_search_t *search; /**< Search (if generated by search). */
    dict_t *dict;	 /**< Dictionary for this DAG. */
    int32 silence;       /**< Silence word ID. */
    int32 frate;         /**< Frame rate. */

    ps_latnode_t *nodes;  /**< List of all nodes. */
    ps_latnode_t *start;  /**< Starting node. */
    ps_latnode_t *end;    /**< Ending node. */

    frame_idx_t n_frames;    /**< Number of frames for this utterance. */
    int16 n_nodes;     /**< Number of nodes in this lattice. */
    int32 final_node_ascr; /**< Acoustic score of implicit link exiting final node. */
    int32 norm;        /**< Normalizer for posterior probabilities. */
    char *hyp_str;     /**< Current hypothesis string. */

    listelem_alloc_t *latnode_alloc;     /**< Node allocator for this DAG. */
    listelem_alloc_t *latlink_alloc;     /**< Link allocator for this DAG. */
    listelem_alloc_t *latlink_list_alloc; /**< List element allocator for this DAG. */

    /* This will probably be replaced with a heap. */
    latlink_list_t *q_head; /**< Queue of links for traversal. */
    latlink_list_t *q_tail; /**< Queue of links for traversal. */
};

/**
 * Links between DAG nodes.
 *
 * A link corresponds to a single hypothesized instance of a word with
 * a given start and end point.

 */
struct ps_latlink_s {
    struct ps_latnode_s *from;	/**< From node */
    struct ps_latnode_s *to;	/**< To node */
    struct ps_latlink_s *best_prev;
    int32 ascr;			/**< Score for from->wid (from->sf to this->ef) */
    int32 path_scr;		/**< Best path score from root of DAG */
    frame_idx_t ef;			/**< Ending frame of this word  */
    int32 alpha;                /**< Forward probability of this link P(w,o_1^{ef}) */
    int32 beta;                 /**< Backward probability of this link P(w|o_{ef+1}^T) */
};

/**
 * DAG nodes.
 *
 * A node corresponds to a number of hypothesized instances of a word
 * which all share the same starting point.
 */
struct ps_latnode_s {
    int32 id;			/**< Unique id for this node */
    int32 wid;			/**< Dictionary word id */
    int32 basewid;		/**< Dictionary base word id */
    /* FIXME: These are (ab)used to store backpointer indices, therefore they MUST be 32 bits. */
    int32 fef;			/**< First end frame */
    int32 lef;			/**< Last end frame */
    frame_idx_t sf;			/**< Start frame */
    int16 reachable;		/**< From \verbatim </s> \endverbatim or \verbatim <s> \endverbatim */
    int32 node_id;		/**< Node from fsg model, used to map lattice back to model */
    union {
        glist_t velist;         /**< List of history entries with different lmstate (tst only) */
	int32 fanin;		/**< Number nodes with links to this node */
	int32 rem_score;	/**< Estimated best score from node.sf to end */
	int32 best_exit;	/**< Best exit score (used for final nodes only) */
    } info;
    latlink_list_t *exits;      /**< Links out of this node */
    latlink_list_t *entries;    /**< Links into this node */

    struct ps_latnode_s *alt;   /**< Node with alternate pronunciation for this word */
    struct ps_latnode_s *next;	/**< Next node in DAG (no ordering implied) */
};

/**
 * Segmentation "iterator" for backpointer table results.
 */
typedef struct dag_seg_s {
    ps_seg_t base;       /**< Base structure. */
    ps_latlink_t **links;   /**< Array of lattice links. */
    int32 norm;     /**< Normalizer for posterior probabilities. */
    int16 n_links;  /**< Number of lattice links. */
    int16 cur;      /**< Current position in bpidx. */
} dag_seg_t;

/**
 * Partial path structure used in N-best (A*) search.
 *
 * Each partial path (latpath_t) is constructed by extending another
 * partial path--parent--by one node.
 */
typedef struct ps_latpath_s {
    ps_latnode_t *node;            /**< Node ending this path. */
    struct ps_latpath_s *parent;   /**< Previous element in this path. */
    struct ps_latpath_s *next;     /**< Pointer to next path in list of paths. */
    int32 score;                  /**< Exact score from start node up to node->sf. */
} ps_latpath_t;

/**
 * A* search structure.
 */
typedef struct ps_astar_s {
    ps_lattice_t *dag;
    ngram_model_t *lmset;
    float32 lwf;

    frame_idx_t sf;
    frame_idx_t ef;
    int32 w1;
    int32 w2;

    int32 n_hyp_tried;
    int32 n_hyp_insert;
    int32 n_hyp_reject;
    int32 insert_depth;
    int32 n_path;

    ps_latpath_t *path_list;
    ps_latpath_t *path_tail;
    ps_latpath_t *top;

    glist_t hyps;	             /**< List of hypothesis strings. */
    listelem_alloc_t *latpath_alloc; /**< Path allocator for N-best search. */
} ps_astar_t;

/**
 * Segmentation "iterator" for A* search results.
 */
typedef struct astar_seg_s {
    ps_seg_t base;
    ps_latnode_t **nodes;
    int n_nodes;
    int cur;
} astar_seg_t;

/**
 * Construct an empty word graph with reference to a search structure.
 */
ps_lattice_t *ps_lattice_init_search(ps_search_t *search, int n_frame);

/**
 * Bypass filler words.
 */
void ps_lattice_bypass_fillers(ps_lattice_t *dag, int32 silpen, int32 fillpen);

/**
 * Remove nodes marked as unreachable.
 */
void ps_lattice_delete_unreachable(ps_lattice_t *dag);

/**
 * Add an edge to the traversal queue.
 */
void ps_lattice_pushq(ps_lattice_t *dag, ps_latlink_t *link);

/**
 * Remove an edge from the traversal queue.
 */
ps_latlink_t *ps_lattice_popq(ps_lattice_t *dag);

/**
 * Clear and reset the traversal queue.
 */
void ps_lattice_delq(ps_lattice_t *dag);

/**
 * Create a new lattice link element.
 */
latlink_list_t *latlink_list_new(ps_lattice_t *dag, ps_latlink_t *link,
                                 latlink_list_t *next);

/**
 * Get hypothesis string after bestpath search.
 */
char const *ps_lattice_hyp(ps_lattice_t *dag, ps_latlink_t *link);

/**
 * Get hypothesis segmentation iterator after bestpath search.
 */
ps_seg_t *ps_lattice_seg_iter(ps_lattice_t *dag, ps_latlink_t *link,
                              float32 lwf);

/**
 * Begin N-Gram based A* search on a word graph.
 *
 * @param sf Starting frame for N-best search.
 * @param ef Ending frame for N-best search, or -1 for last frame.
 * @param w1 First context word, or -1 for none.
 * @param w2 Second context word, or -1 for none.
 * @return 0 for success, <0 on error.
 */
ps_astar_t *ps_astar_start(ps_lattice_t *dag,
                           ngram_model_t *lmset,
                           float32 lwf,
                           int sf, int ef,
                           int w1, int w2);

/**
 * Find next best hypothesis of A* on a word graph.
 *
 * @return a complete path, or NULL if no more hypotheses exist.
 */
ps_latpath_t *ps_astar_next(ps_astar_t *nbest);

/**
 * Finish N-best search, releasing resources associated with it.
 */
void ps_astar_finish(ps_astar_t *nbest);

/**
 * Get hypothesis string from A* search.
 */
char const *ps_astar_hyp(ps_astar_t *nbest, ps_latpath_t *path);

/**
 * Get hypothesis segmentation from A* search.
 */
ps_seg_t *ps_astar_seg_iter(ps_astar_t *astar, ps_latpath_t *path, float32 lwf);


#endif /* __PS_LATTICE_INTERNAL_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2009 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file ps_mllr.c Model-space linear transforms for speaker adaptation
 */

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>

/* Local headers. */
#include "acmod.h"

ps_mllr_t *
ps_mllr_read(char const *regmatfile)
{
    ps_mllr_t *mllr;
    FILE *fp;
    int n, i, m, j, k;

    mllr = ckd_calloc(1, sizeof(*mllr));
    mllr->refcnt = 1;

    if ((fp = fopen(regmatfile, "r")) == NULL) {
        E_ERROR_SYSTEM("Failed to open MLLR file '%s' for reading", regmatfile);
        goto error_out;
    }
    else
        E_INFO("Reading MLLR transformation file '%s'\n", regmatfile);

    if ((fscanf(fp, "%d", &n) != 1) || (n < 1)) {
        E_ERROR("Failed to read number of MLLR classes\n");
        goto error_out;
    }
    mllr->n_class = n;

    if ((fscanf(fp, "%d", &n) != 1)) {
        E_ERROR("Failed to read number of feature streams\n");
        goto error_out;
    }
    mllr->n_feat = n;
    mllr->veclen = ckd_calloc(mllr->n_feat, sizeof(*mllr->veclen));

    mllr->A = (float32 ****) ckd_calloc(mllr->n_feat, sizeof(float32 **));
    mllr->b = (float32 ***) ckd_calloc(mllr->n_feat, sizeof(float32 *));
    mllr->h = (float32 ***) ckd_calloc(mllr->n_feat, sizeof(float32 *));

    for (i = 0; i < mllr->n_feat; ++i) {
        if (fscanf(fp, "%d", &n) != 1) {
            E_ERROR("Failed to read stream length for feature %d\n", i);
            goto error_out;
        }
        mllr->veclen[i] = n;
        mllr->A[i] =
            (float32 ***) ckd_calloc_3d(mllr->n_class, mllr->veclen[i],
                                        mllr->veclen[i], sizeof(float32));
        mllr->b[i] =
            (float32 **) ckd_calloc_2d(mllr->n_class, mllr->veclen[i],
                                       sizeof(float32));
        mllr->h[i] =
            (float32 **) ckd_calloc_2d(mllr->n_class, mllr->veclen[i],
                                       sizeof(float32));

        for (m = 0; m < mllr->n_class; ++m) {
            for (j = 0; j < mllr->veclen[i]; ++j) {
                for (k = 0; k < mllr->veclen[i]; ++k) {
                    if (fscanf(fp, "%f ", &mllr->A[i][m][j][k]) != 1) {
                        E_ERROR("Failed reading MLLR rotation (%d,%d,%d,%d)\n",
                                i, m, j, k);
                        goto error_out;
                    }
                }
            }
            for (j = 0; j < mllr->veclen[i]; ++j) {
                if (fscanf(fp, "%f ", &mllr->b[i][m][j]) != 1) {
                    E_ERROR("Failed reading MLLR bias (%d,%d,%d)\n",
                            i, m, j);
                    goto error_out;
                }
            }
            for (j = 0; j < mllr->veclen[i]; ++j) {
                if (fscanf(fp, "%f ", &mllr->h[i][m][j]) != 1) {
                    E_ERROR("Failed reading MLLR variance scale (%d,%d,%d)\n",
                            i, m, j);
                    goto error_out;
                }
            }
        }
    }
    fclose(fp);
    return mllr;

error_out:
    if (fp)
        fclose(fp);
    ps_mllr_free(mllr);
    return NULL;
}

ps_mllr_t *
ps_mllr_retain(ps_mllr_t *mllr)
{
    ++mllr->refcnt;
    return mllr;
}

int
ps_mllr_free(ps_mllr_t *mllr)
{
    int i;

    if (mllr == NULL)
        return 0;
    if (--mllr->refcnt > 0)
        return mllr->refcnt;

    for (i = 0; i < mllr->n_feat; ++i) {
        if (mllr->A)
            ckd_free_3d(mllr->A[i]);
        if (mllr->b)
            ckd_free_2d(mllr->b[i]);
        if (mllr->h)
            ckd_free_2d(mllr->h[i]);
    }
    ckd_free(mllr->veclen);
    ckd_free(mllr->A);
    ckd_free(mllr->b);
    ckd_free(mllr->h);
    ckd_free(mllr);

    return 0;
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#if defined(__ADSPBLACKFIN__)
#elif !defined(_WIN32_WCE)
#include <sys/types.h>
#endif

/* SphinxBase headers */
#include <sphinx_config.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/fixpoint.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/bio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/prim_type.h>

/* Local headers */
#include "tied_mgau_common.h"
#include "ptm_mgau.h"

static ps_mgaufuncs_t ptm_mgau_funcs = {
    "ptm",
    ptm_mgau_frame_eval,      /* frame_eval */
    ptm_mgau_mllr_transform,  /* transform */
    ptm_mgau_free             /* free */
};

#define COMPUTE_GMM_MAP(_idx)                           \
    diff[_idx] = obs[_idx] - mean[_idx];                \
    sqdiff[_idx] = MFCCMUL(diff[_idx], diff[_idx]);     \
    compl[_idx] = MFCCMUL(sqdiff[_idx], var[_idx]);
#define COMPUTE_GMM_REDUCE(_idx)                \
    d = GMMSUB(d, compl[_idx]);

static void
insertion_sort_topn(ptm_topn_t *topn, int i, int32 d)
{
    ptm_topn_t vtmp;
    int j;

    topn[i].score = d;
    if (i == 0)
        return;
    vtmp = topn[i];
    for (j = i - 1; j >= 0 && d > topn[j].score; j--) {
        topn[j + 1] = topn[j];
    }
    topn[j + 1] = vtmp;
}

static int
eval_topn_one(ptm_mgau_t *s, int cb, int feat, mfcc_t *z) //There are 2 eval topn functions defined.First one is this andisrenamed eval_topn_one
{
    ptm_topn_t *topn;
    int i, ceplen;

    topn = s->f->topn[cb][feat];
    ceplen = s->g->featlen[feat];

    for (i = 0; i < s->max_topn; i++) {
        mfcc_t *mean, diff[4], sqdiff[4], compl[4]; /* diff, diff^2, component likelihood */
        mfcc_t *var, d;
        mfcc_t *obs;
        int32 cw, j;

        cw = topn[i].cw;
        mean = s->g->mean[cb][feat][0] + cw * ceplen;
        var = s->g->var[cb][feat][0] + cw * ceplen;
        d = s->g->det[cb][feat][cw];
        obs = z;
        for (j = 0; j < ceplen % 4; ++j) {
            diff[0] = *obs++ - *mean++;
            sqdiff[0] = MFCCMUL(diff[0], diff[0]);
            compl[0] = MFCCMUL(sqdiff[0], *var);
            d = GMMSUB(d, compl[0]);
            ++var;
        }
        /* We could vectorize this but it's unlikely to make much
         * difference as the outer loop here isn't very big. */
        for (;j < ceplen; j += 4) {
            COMPUTE_GMM_MAP(0);
            COMPUTE_GMM_MAP(1);
            COMPUTE_GMM_MAP(2);
            COMPUTE_GMM_MAP(3);
            COMPUTE_GMM_REDUCE(0);
            COMPUTE_GMM_REDUCE(1);
            COMPUTE_GMM_REDUCE(2);
            COMPUTE_GMM_REDUCE(3);
            var += 4;
            obs += 4;
            mean += 4;
        }
        insertion_sort_topn(topn, i, (int32)d);
    }

    return topn[0].score;
}

/* This looks bad, but it actually isn't.  Less than 1% of eval_cb's
 * time is spent doing this. */
static void
insertion_sort_cb(ptm_topn_t **cur, ptm_topn_t *worst, ptm_topn_t *best,
                  int cw, int32 intd)
{
    for (*cur = worst - 1; *cur >= best && intd >= (*cur)->score; --*cur)
        memcpy(*cur + 1, *cur, sizeof(**cur));
    ++*cur;
    (*cur)->cw = cw;
    (*cur)->score = intd;
}

static int
eval_cb_one(ptm_mgau_t *s, int cb, int feat, mfcc_t *z) //#######################################################################
{
    ptm_topn_t *worst, *best, *topn;
    mfcc_t *mean;
    mfcc_t *var, *det, *detP, *detE;
    int32 i, ceplen;

    best = topn = s->f->topn[cb][feat];
    worst = topn + (s->max_topn - 1);
    mean = s->g->mean[cb][feat][0];
    var = s->g->var[cb][feat][0];
    det = s->g->det[cb][feat];
    detE = det + s->g->n_density;
    ceplen = s->g->featlen[feat];

    for (detP = det; detP < detE; ++detP) {
        mfcc_t diff[4], sqdiff[4], compl[4]; /* diff, diff^2, component likelihood */
        mfcc_t d, thresh;
        mfcc_t *obs;
        ptm_topn_t *cur;
        int32 cw, j;

        d = *detP;
        thresh = (mfcc_t) worst->score; /* Avoid int-to-float conversions */
        obs = z;
        cw = detP - det;

        /* Unroll the loop starting with the first dimension(s).  In
         * theory this might be a bit faster if this Gaussian gets
         * "knocked out" by C0. In practice not. */
        for (j = 0; (j < ceplen % 4) && (d >= thresh); ++j) {
            diff[0] = *obs++ - *mean++;
            sqdiff[0] = MFCCMUL(diff[0], diff[0]);
            compl[0] = MFCCMUL(sqdiff[0], *var++);
            d = GMMSUB(d, compl[0]);
        }
        /* Now do 4 dimensions at a time.  You'd think that GCC would
         * vectorize this?  Apparently not.  And it's right, because
         * that won't make this any faster, at least on x86-64. */
        for (; j < ceplen && d >= thresh; j += 4) {
            COMPUTE_GMM_MAP(0);
            COMPUTE_GMM_MAP(1);
            COMPUTE_GMM_MAP(2);
            COMPUTE_GMM_MAP(3);
            COMPUTE_GMM_REDUCE(0);
            COMPUTE_GMM_REDUCE(1);
            COMPUTE_GMM_REDUCE(2);
            COMPUTE_GMM_REDUCE(3);
            var += 4;
            obs += 4;
            mean += 4;
        }
        if (j < ceplen) {
            /* terminated early, so not in topn */
            mean += (ceplen - j);
            var += (ceplen - j);
            continue;
        }
        if (d < thresh)
            continue;
        for (i = 0; i < s->max_topn; i++) {
            /* already there, so don't need to insert */
            if (topn[i].cw == cw)
                break;
        }
        if (i < s->max_topn)
            continue;       /* already there.  Don't insert */
        insertion_sort_cb(&cur, worst, best, cw, (int32)d);
    }

    return best->score;
}

/**
 * Compute top-N densities for active codebooks (and prune)
 */
static int
ptm_mgau_codebook_eval(ptm_mgau_t *s, mfcc_t **z, int frame)
{
    int i, j;

    /* First evaluate top-N from previous frame. */
    for (i = 0; i < s->g->n_mgau; ++i)
        for (j = 0; j < s->g->n_feat; ++j)
            eval_topn_one(s, i, j, z[j]);

    /* If frame downsampling is in effect, possibly do nothing else. */
    if (frame % s->ds_ratio)
        return 0;

    /* Evaluate remaining codebooks. */
    for (i = 0; i < s->g->n_mgau; ++i) {
        if (bitvec_is_clear(s->f->mgau_active, i))
            continue;
        for (j = 0; j < s->g->n_feat; ++j) {
            eval_cb_one(s, i, j, z[j]);
        }
    }

    /* Normalize densities to produce "posterior probabilities",
     * i.e. things with a reasonable dynamic range, then scale and
     * clamp them to the acceptable range.  This is actually done
     * solely to ensure that we can use fast_logmath_add().  Note that
     * unless we share the same normalizer across all codebooks for
     * each feature stream we get defective scores (that's why these
     * loops are inside out - doing it per-feature should give us
     * greater precision). */
    for (j = 0; j < s->g->n_feat; ++j) {
        int32 norm = 0x7fffffff;
        for (i = 0; i < s->g->n_mgau; ++i) {
            if (bitvec_is_clear(s->f->mgau_active, i))
                continue;
            if (norm > s->f->topn[i][j][0].score >> SENSCR_SHIFT)
                norm = s->f->topn[i][j][0].score >> SENSCR_SHIFT;
        }
        assert(norm != 0x7fffffff);
        for (i = 0; i < s->g->n_mgau; ++i) {
            int32 k;
            if (bitvec_is_clear(s->f->mgau_active, i))
                continue;
            for (k = 0; k < s->max_topn; ++k) {
                s->f->topn[i][j][k].score >>= SENSCR_SHIFT;
                s->f->topn[i][j][k].score -= norm;
                s->f->topn[i][j][k].score = -s->f->topn[i][j][k].score;
                if (s->f->topn[i][j][k].score > MAX_NEG_ASCR) 
                    s->f->topn[i][j][k].score = MAX_NEG_ASCR;
            }
        }
    }

    return 0;
}

static int
ptm_mgau_calc_cb_active(ptm_mgau_t *s, uint8 *senone_active,
                        int32 n_senone_active, int compallsen)
{
    int i, lastsen;

    if (compallsen) {
        bitvec_set_all(s->f->mgau_active, s->g->n_mgau);
        return 0;
    }
    bitvec_clear_all(s->f->mgau_active, s->g->n_mgau);
    for (lastsen = i = 0; i < n_senone_active; ++i) {
        int sen = senone_active[i] + lastsen;
        int cb = s->sen2cb[sen];
        bitvec_set(s->f->mgau_active, cb);
        lastsen = sen;
    }
    E_DEBUG(1, ("Active codebooks:"));
    for (i = 0; i < s->g->n_mgau; ++i) {
        if (bitvec_is_clear(s->f->mgau_active, i))
            continue;
        E_DEBUGCONT(1, (" %d", i));
    }
    E_DEBUGCONT(1, ("\n"));
    return 0;
}

/**
 * Compute senone scores from top-N densities for active codebooks.
 */
static int
ptm_mgau_senone_eval(ptm_mgau_t *s, int16 *senone_scores,
                     uint8 *senone_active, int32 n_senone_active,
                     int compall)
{
    int i, lastsen, bestscore;

    memset(senone_scores, 0, s->n_sen * sizeof(*senone_scores));
    /* FIXME: This is the non-cache-efficient way to do this.  We want
     * to evaluate one codeword at a time but this requires us to have
     * a reverse codebook to senone mapping, which we don't have
     * (yet), since different codebooks have different top-N
     * codewords. */
    if (compall)
        n_senone_active = s->n_sen;
    bestscore = 0x7fffffff;
    for (lastsen = i = 0; i < n_senone_active; ++i) {
        int sen, f, cb;
        int ascore;

        if (compall)
            sen = i;
        else
            sen = senone_active[i] + lastsen;
        lastsen = sen;
        cb = s->sen2cb[sen];

        if (bitvec_is_clear(s->f->mgau_active, cb)) {
            int j;
            /* Because senone_active is deltas we can't really "knock
             * out" senones from pruned codebooks, and in any case,
             * it wouldn't make any difference to the search code,
             * which doesn't expect senone_active to change. */
            for (f = 0; f < s->g->n_feat; ++f) {
                for (j = 0; j < s->max_topn; ++j) {
                    s->f->topn[cb][f][j].score = MAX_NEG_ASCR;
                }
            }
        }
        /* For each feature, log-sum codeword scores + mixw to get
         * feature density, then sum (multiply) to get ascore */
        ascore = 0;
        for (f = 0; f < s->g->n_feat; ++f) {
            ptm_topn_t *topn;
            int j, fden = 0;
            topn = s->f->topn[cb][f];
            for (j = 0; j < s->max_topn; ++j) {
                int mixw;
                /* Find mixture weight for this codeword. */
                if (s->mixw_cb) {
                    int dcw = s->mixw[f][topn[j].cw][sen/2];
                    dcw = (dcw & 1) ? dcw >> 4 : dcw & 0x0f;
                    mixw = s->mixw_cb[dcw];
                }
                else {
                    mixw = s->mixw[f][topn[j].cw][sen];
                }
                if (j == 0)
                    fden = mixw + topn[j].score;
                else
                    fden = fast_logmath_add(s->lmath_8b, fden,
                                       mixw + topn[j].score);
                E_DEBUG(3, ("fden[%d][%d] l+= %d + %d = %d\n",
                            sen, f, mixw, topn[j].score, fden));
            }
            ascore += fden;
        }
        if (ascore < bestscore) bestscore = ascore;
        senone_scores[sen] = ascore;
    }
    /* Normalize the scores again (finishing the job we started above
     * in ptm_mgau_codebook_eval...) */
    for (i = 0; i < s->n_sen; ++i) {
        senone_scores[i] -= bestscore;
    }

    return 0;
}

/**
 * Compute senone scores for the active senones.
 */
int32
ptm_mgau_frame_eval(ps_mgau_t *ps,
                    int16 *senone_scores,
                    uint8 *senone_active,
                    int32 n_senone_active,
                    mfcc_t ** featbuf, int32 frame,
                    int32 compallsen)
{
    ptm_mgau_t *s = (ptm_mgau_t *)ps;
    int fast_eval_idx;

    /* Find the appropriate frame in the rotating history buffer
     * corresponding to the requested input frame.  No bounds checking
     * is done here, which just means you'll get semi-random crap if
     * you request a frame in the future or one that's too far in the
     * past.  Since the history buffer is just used for fast match
     * that might not be fatal. */
    fast_eval_idx = frame % s->n_fast_hist;
    s->f = s->hist + fast_eval_idx;
    /* Compute the top-N codewords for every codebook, unless this
     * is a past frame, in which case we already have them (we
     * hope!) */
    if (frame >= ps_mgau_base(ps)->frame_idx) {
        ptm_fast_eval_t *lastf;
        /* Get the previous frame's top-N information (on the
         * first frame of the input this is just all WORST_DIST,
         * no harm in that) */
        if (fast_eval_idx == 0)
            lastf = s->hist + s->n_fast_hist - 1;
        else
            lastf = s->hist + fast_eval_idx - 1;
        /* Copy in initial top-N info */
        memcpy(s->f->topn[0][0], lastf->topn[0][0],
               s->g->n_mgau * s->g->n_feat * s->max_topn * sizeof(ptm_topn_t));
        /* Generate initial active codebook list (this might not be
         * necessary) */
        ptm_mgau_calc_cb_active(s, senone_active, n_senone_active, compallsen);
        /* Now evaluate top-N, prune, and evaluate remaining codebooks. */
        ptm_mgau_codebook_eval(s, featbuf, frame);
    }
    /* Evaluate intersection of active senones and active codebooks. */
    ptm_mgau_senone_eval(s, senone_scores, senone_active,
                         n_senone_active, compallsen);

    return 0;
}

static int32
read_sendump_one(ptm_mgau_t *s, bin_mdef_t *mdef, char const *file)
{
    FILE *fp;
    char line[1000];
    int32 i, n, r, c;
    int32 do_swap, do_mmap;
    size_t offset;
    int n_clust = 0;
    int n_feat = s->g->n_feat;
    int n_density = s->g->n_density;
    int n_sen = bin_mdef_n_sen(mdef);
    int n_bits = 8;

    s->n_sen = n_sen; /* FIXME: Should have been done earlier */
    do_mmap = cmd_ln_boolean_r(s->config, "-mmap");

    if ((fp = fopen(file, "rb")) == NULL)
        return -1;

    E_INFO("Loading senones from dump file %s\n", file);
    /* Read title size, title */
    if (fread(&n, sizeof(int32), 1, fp) != 1) {
        E_ERROR_SYSTEM("Failed to read title size from %s", file);
        goto error_out;
    }
    /* This is extremely bogus */
    do_swap = 0;
    if (n < 1 || n > 999) {
        SWAP_INT32(&n);
        if (n < 1 || n > 999) {
            E_ERROR("Title length %x in dump file %s out of range\n", n, file);
            goto error_out;
        }
        do_swap = 1;
    }
    if (fread(line, sizeof(char), n, fp) != n) {
        E_ERROR_SYSTEM("Cannot read title");
        goto error_out;
    }
    if (line[n - 1] != '\0') {
        E_ERROR("Bad title in dump file\n");
        goto error_out;
    }
    E_INFO("%s\n", line);

    /* Read header size, header */
    if (fread(&n, sizeof(n), 1, fp) != 1) {
        E_ERROR_SYSTEM("Failed to read header size from %s", file);
        goto error_out;
    }
    if (do_swap) SWAP_INT32(&n);
    if (fread(line, sizeof(char), n, fp) != n) {
        E_ERROR_SYSTEM("Cannot read header");
        goto error_out;
    }
    if (line[n - 1] != '\0') {
        E_ERROR("Bad header in dump file\n");
        goto error_out;
    }

    /* Read other header strings until string length = 0 */
    for (;;) {
        if (fread(&n, sizeof(n), 1, fp) != 1) {
            E_ERROR_SYSTEM("Failed to read header string size from %s", file);
            goto error_out;
        }
        if (do_swap) SWAP_INT32(&n);
        if (n == 0)
            break;
        if (fread(line, sizeof(char), n, fp) != n) {
            E_ERROR_SYSTEM("Cannot read header");
            goto error_out;
        }
        /* Look for a cluster count, if present */
        if (!strncmp(line, "feature_count ", strlen("feature_count "))) {
            n_feat = atoi(line + strlen("feature_count "));
        }
        if (!strncmp(line, "mixture_count ", strlen("mixture_count "))) {
            n_density = atoi(line + strlen("mixture_count "));
        }
        if (!strncmp(line, "model_count ", strlen("model_count "))) {
            n_sen = atoi(line + strlen("model_count "));
        }
        if (!strncmp(line, "cluster_count ", strlen("cluster_count "))) {
            n_clust = atoi(line + strlen("cluster_count "));
        }
        if (!strncmp(line, "cluster_bits ", strlen("cluster_bits "))) {
            n_bits = atoi(line + strlen("cluster_bits "));
        }
    }

    /* Defaults for #rows, #columns in mixw array. */
    c = n_sen;
    r = n_density;
    if (n_clust == 0) {
        /* Older mixw files have them here, and they might be padded. */
        if (fread(&r, sizeof(r), 1, fp) != 1) {
            E_ERROR_SYSTEM("Cannot read #rows");
            goto error_out;
        }
        if (do_swap) SWAP_INT32(&r);
        if (fread(&c, sizeof(c), 1, fp) != 1) {
            E_ERROR_SYSTEM("Cannot read #columns");
            goto error_out;
        }
        if (do_swap) SWAP_INT32(&c);
        E_INFO("Rows: %d, Columns: %d\n", r, c);
    }

    if (n_feat != s->g->n_feat) {
        E_ERROR("Number of feature streams mismatch: %d != %d\n",
                n_feat, s->g->n_feat);
        goto error_out;
    }
    if (n_density != s->g->n_density) {
        E_ERROR("Number of densities mismatch: %d != %d\n",
                n_density, s->g->n_density);
        goto error_out;
    }
    if (n_sen != s->n_sen) {
        E_ERROR("Number of senones mismatch: %d != %d\n",
                n_sen, s->n_sen);
        goto error_out;
    }

    if (!((n_clust == 0) || (n_clust == 15) || (n_clust == 16))) {
        E_ERROR("Cluster count must be 0, 15, or 16\n");
        goto error_out;
    }
    if (n_clust == 15)
        ++n_clust;

    if (!((n_bits == 8) || (n_bits == 4))) {
        E_ERROR("Cluster count must be 4 or 8\n");
        goto error_out;
    }

    if (do_mmap) {
            E_INFO("Using memory-mapped I/O for senones\n");
    }
    offset = ftell(fp);

    /* Allocate memory for pdfs (or memory map them) */
    if (do_mmap) {
        s->sendump_mmap = mmio_file_read(file);
        /* Get cluster codebook if any. */
        if (n_clust) {
            s->mixw_cb = ((uint8 *) mmio_file_ptr(s->sendump_mmap)) + offset;
            offset += n_clust;
        }
    }
    else {
        /* Get cluster codebook if any. */
        if (n_clust) {
            s->mixw_cb = ckd_calloc(1, n_clust);
            if (fread(s->mixw_cb, 1, n_clust, fp) != (size_t) n_clust) {
                E_ERROR("Failed to read %d bytes from sendump\n", n_clust);
                goto error_out;
            }
        }
    }

    /* Set up pointers, or read, or whatever */
    if (s->sendump_mmap) {
        s->mixw = ckd_calloc_2d(n_feat, n_density, sizeof(*s->mixw));
        for (n = 0; n < n_feat; n++) {
            int step = c;
            if (n_bits == 4)
                step = (step + 1) / 2;
            for (i = 0; i < r; i++) {
                s->mixw[n][i] = ((uint8 *) mmio_file_ptr(s->sendump_mmap)) + offset;
                offset += step;
            }
        }
    }
    else {
        s->mixw = ckd_calloc_3d(n_feat, n_density, n_sen, sizeof(***s->mixw));
        /* Read pdf values and ids */
        for (n = 0; n < n_feat; n++) {
            int step = c;
            if (n_bits == 4)
                step = (step + 1) / 2;
            for (i = 0; i < r; i++) {
                if (fread(s->mixw[n][i], sizeof(***s->mixw), step, fp)
                    != (size_t) step) {
                    E_ERROR("Failed to read %d bytes from sendump\n", step);
                    goto error_out;
                }
            }
        }
    }

    fclose(fp);
    return 0;
error_out:
    fclose(fp);
    return -1;
}

static int32
read_mixw_one(ptm_mgau_t * s, char const *file_name, double SmoothMin)
{
    char **argname, **argval;
    char eofchk;
    FILE *fp;
    int32 byteswap, chksum_present;
    uint32 chksum;
    float32 *pdf;
    int32 i, f, c, n;
    int32 n_sen;
    int32 n_feat;
    int32 n_comp;
    int32 n_err;

    E_INFO("Reading mixture weights file '%s'\n", file_name);

    if ((fp = fopen(file_name, "rb")) == NULL)
        E_FATAL_SYSTEM("Failed to open mixture file '%s' for reading", file_name);

    /* Read header, including argument-value info and 32-bit byteorder magic */
    if (bio_readhdr(fp, &argname, &argval, &byteswap) < 0)
        E_FATAL("Failed to read header from '%s'\n", file_name);

    /* Parse argument-value list */
    chksum_present = 0;
    for (i = 0; argname[i]; i++) {
        if (strcmp(argname[i], "version") == 0) {
            if (strcmp(argval[i], MGAU_MIXW_VERSION) != 0)
                E_WARN("Version mismatch(%s): %s, expecting %s\n",
                       file_name, argval[i], MGAU_MIXW_VERSION);
        }
        else if (strcmp(argname[i], "chksum0") == 0) {
            chksum_present = 1; /* Ignore the associated value */
        }
    }
    bio_hdrarg_free(argname, argval);
    argname = argval = NULL;

    chksum = 0;

    /* Read #senones, #features, #codewords, arraysize */
    if ((bio_fread(&n_sen, sizeof(int32), 1, fp, byteswap, &chksum) != 1)
        || (bio_fread(&n_feat, sizeof(int32), 1, fp, byteswap, &chksum) !=
            1)
        || (bio_fread(&n_comp, sizeof(int32), 1, fp, byteswap, &chksum) !=
            1)
        || (bio_fread(&n, sizeof(int32), 1, fp, byteswap, &chksum) != 1)) {
        E_FATAL("bio_fread(%s) (arraysize) failed\n", file_name);
    }
    if (n_feat != s->g->n_feat)
        E_FATAL("#Features streams(%d) != %d\n", n_feat, s->g->n_feat);
    if (n != n_sen * n_feat * n_comp) {
        E_FATAL
            ("%s: #float32s(%d) doesn't match header dimensions: %d x %d x %d\n",
             file_name, i, n_sen, n_feat, n_comp);
    }

    /* n_sen = number of mixture weights per codeword, which is
     * fixed at the number of senones since we have only one codebook.
     */
    s->n_sen = n_sen;

    /* Quantized mixture weight arrays. */
    s->mixw = ckd_calloc_3d(s->g->n_feat, s->g->n_density,
                            n_sen, sizeof(***s->mixw));

    /* Temporary structure to read in floats before conversion to (int32) logs3 */
    pdf = (float32 *) ckd_calloc(n_comp, sizeof(float32));

    /* Read senone probs data, normalize, floor, convert to logs3, truncate to 8 bits */
    n_err = 0;
    for (i = 0; i < n_sen; i++) {
        for (f = 0; f < n_feat; f++) {
            if (bio_fread((void *) pdf, sizeof(float32),
                          n_comp, fp, byteswap, &chksum) != n_comp) {
                E_FATAL("bio_fread(%s) (arraydata) failed\n", file_name);
            }

            /* Normalize and floor */
            if (vector_sum_norm(pdf, n_comp) <= 0.0)
                n_err++;
            vector_floor(pdf, n_comp, SmoothMin);
            vector_sum_norm(pdf, n_comp);

            /* Convert to LOG, quantize, and transpose */
            for (c = 0; c < n_comp; c++) {
                int32 qscr;

                qscr = -logmath_log(s->lmath_8b, pdf[c]);
                if ((qscr > MAX_NEG_MIXW) || (qscr < 0))
                    qscr = MAX_NEG_MIXW;
                s->mixw[f][c][i] = qscr;
            }
        }
    }
    if (n_err > 0)
        E_WARN("Weight normalization failed for %d mixture weights components\n", n_err);

    ckd_free(pdf);

    if (chksum_present)
        bio_verify_chksum(fp, byteswap, chksum);

    if (fread(&eofchk, 1, 1, fp) == 1)
        E_FATAL("More data than expected in %s\n", file_name);

    fclose(fp);

    E_INFO("Read %d x %d x %d mixture weights\n", n_sen, n_feat, n_comp);
    return n_sen;
}

ps_mgau_t *
ptm_mgau_init(acmod_t *acmod, bin_mdef_t *mdef)
{
    ptm_mgau_t *s;
    ps_mgau_t *ps;
    char const *sendump_path;
    int i;

    s = ckd_calloc(1, sizeof(*s));
    s->config = acmod->config;

    s->lmath = logmath_retain(acmod->lmath);
    /* Log-add table. */
    s->lmath_8b = logmath_init(logmath_get_base(acmod->lmath), SENSCR_SHIFT, TRUE);
    if (s->lmath_8b == NULL)
        goto error_out;
    /* Ensure that it is only 8 bits wide so that fast_logmath_add() works. */
    if (logmath_get_width(s->lmath_8b) != 1) {
        E_ERROR("Log base %f is too small to represent add table in 8 bits\n",
                logmath_get_base(s->lmath_8b));
        goto error_out;
    }

    /* Read means and variances. */
    if ((s->g = gauden_init(cmd_ln_str_r(s->config, "-mean"),
                            cmd_ln_str_r(s->config, "-var"),
                            cmd_ln_float32_r(s->config, "-varfloor"),
                            s->lmath)) == NULL)
        goto error_out;
    /* We only support 256 codebooks or less (like 640k or 2GB, this
     * should be enough for anyone) */
    if (s->g->n_mgau > 256) {
        E_INFO("Number of codebooks exceeds 256: %d\n", s->g->n_mgau);
        goto error_out;
    }
    if (s->g->n_mgau != bin_mdef_n_ciphone(mdef)) {
        E_INFO("Number of codebooks doesn't match number of ciphones, doesn't look like PTM: %d != %d\n", s->g->n_mgau, bin_mdef_n_ciphone(mdef));
        goto error_out;
    }
    /* Verify n_feat and veclen, against acmod. */
    if (s->g->n_feat != feat_dimension1(acmod->fcb)) {
        E_ERROR("Number of streams does not match: %d != %d\n",
                s->g->n_feat, feat_dimension1(acmod->fcb));
        goto error_out;
    }
    for (i = 0; i < s->g->n_feat; ++i) {
        if (s->g->featlen[i] != feat_dimension2(acmod->fcb, i)) {
            E_ERROR("Dimension of stream %d does not match: %d != %d\n",
                    s->g->featlen[i], feat_dimension2(acmod->fcb, i));
            goto error_out;
        }
    }
    /* Read mixture weights. */
    if ((sendump_path = cmd_ln_str_r(s->config, "-sendump"))) {
        if (read_sendump_one(s, acmod->mdef, sendump_path) < 0) {
            goto error_out;
        }
    }
    else {
        if (read_mixw_one(s, cmd_ln_str_r(s->config, "-mixw"),
                      cmd_ln_float32_r(s->config, "-mixwfloor")) < 0) {
            goto error_out;
        }
    }
    s->ds_ratio = cmd_ln_int32_r(s->config, "-ds");
    s->max_topn = cmd_ln_int32_r(s->config, "-topn");
    E_INFO("Maximum top-N: %d\n", s->max_topn);

    /* Assume mapping of senones to their base phones, though this
     * will become more flexible in the future. */
    s->sen2cb = ckd_calloc(s->n_sen, sizeof(*s->sen2cb));
    for (i = 0; i < s->n_sen; ++i)
        s->sen2cb[i] = bin_mdef_sen2cimap(acmod->mdef, i);

    /* Allocate fast-match history buffers.  We need enough for the
     * phoneme lookahead window, plus the current frame, plus one for
     * good measure? (FIXME: I don't remember why) */
    s->n_fast_hist = cmd_ln_int32_r(s->config, "-pl_window") + 2;
    s->hist = ckd_calloc(s->n_fast_hist, sizeof(*s->hist));
    /* s->f will be a rotating pointer into s->hist. */
    s->f = s->hist;
    for (i = 0; i < s->n_fast_hist; ++i) {
        int j, k, m;
        /* Top-N codewords for every codebook and feature. */
        s->hist[i].topn = ckd_calloc_3d(s->g->n_mgau, s->g->n_feat,
                                        s->max_topn, sizeof(ptm_topn_t));
        /* Initialize them to sane (yet arbitrary) defaults. */
        for (j = 0; j < s->g->n_mgau; ++j) {
            for (k = 0; k < s->g->n_feat; ++k) {
                for (m = 0; m < s->max_topn; ++m) {
                    s->hist[i].topn[j][k][m].cw = m;
                    s->hist[i].topn[j][k][m].score = WORST_DIST;
                }
            }
        }
        /* Active codebook mapping (just codebook, not features,
           at least not yet) */
        s->hist[i].mgau_active = bitvec_alloc(s->g->n_mgau);
        /* Start with them all on, prune them later. */
        bitvec_set_all(s->hist[i].mgau_active, s->g->n_mgau);
    }

    ps = (ps_mgau_t *)s;
    ps->vt = &ptm_mgau_funcs;
    return ps;
error_out:
    ptm_mgau_free(ps_mgau_base(s));
    return NULL;
}

int
ptm_mgau_mllr_transform(ps_mgau_t *ps,
                            ps_mllr_t *mllr)
{
    ptm_mgau_t *s = (ptm_mgau_t *)ps;
    return gauden_mllr_transform(s->g, mllr, s->config);
}

void
ptm_mgau_free(ps_mgau_t *ps)
{
    ptm_mgau_t *s = (ptm_mgau_t *)ps;

    logmath_free(s->lmath);
    logmath_free(s->lmath_8b);
    if (s->sendump_mmap) {
        ckd_free_2d(s->mixw); 
        mmio_file_unmap(s->sendump_mmap);
    }
    else {
        ckd_free_3d(s->mixw);
    }
    ckd_free(s->sen2cb);
    gauden_free(s->g);
    ckd_free(s);
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/**
 * @file ptm_mgau.h Fast phonetically-tied mixture evaluation.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __PTM_MGAU_H__
#define __PTM_MGAU_H__

/* SphinxBase headesr. */
#include <sphinxbase/fe.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/mmio.h>

/* Local headers. */
#include "acmod.h"
#include "hmm.h"
#include "bin_mdef.h"
#include "ms_gauden.h"

typedef struct ptm_mgau_s ptm_mgau_t;

typedef struct ptm_topn_s {
    int32 cw;    /**< Codeword index. */
    int32 score; /**< Score. */
} ptm_topn_t;

typedef struct ptm_fast_eval_s {
    ptm_topn_t ***topn;     /**< Top-N for each codebook (mgau x feature x topn) */
    bitvec_t *mgau_active; /**< Set of active codebooks */
} ptm_fast_eval_t;

struct ptm_mgau_s {
    ps_mgau_t base;     /**< base structure. */
    cmd_ln_t *config;   /**< Configuration parameters */
    gauden_t *g;        /**< Set of Gaussians. */
    int32 n_sen;       /**< Number of senones. */
    uint8 *sen2cb;     /**< Senone to codebook mapping. */
    uint8 ***mixw;     /**< Mixture weight distributions by feature, codeword, senone */
    mmio_file_t *sendump_mmap;/* Memory map for mixw (or NULL if not mmap) */
    uint8 *mixw_cb;    /* Mixture weight codebook, if any (assume it contains 16 values) */
    int16 max_topn;
    int16 ds_ratio;

    ptm_fast_eval_t *hist;   /**< Fast evaluation info for past frames. */
    ptm_fast_eval_t *f;      /**< Fast eval info for current frame. */
    int n_fast_hist;         /**< Number of past frames tracked. */

    /* Log-add table for compressed values. */
    logmath_t *lmath_8b;
    /* Log-add object for reloading means/variances. */
    logmath_t *lmath;
};

ps_mgau_t *ptm_mgau_init(acmod_t *acmod, bin_mdef_t *mdef);
void ptm_mgau_free(ps_mgau_t *s);
int ptm_mgau_frame_eval(ps_mgau_t *s,
                        int16 *senone_scores,
                        uint8 *senone_active,
                        int32 n_senone_active,
                        mfcc_t **featbuf,
                        int32 frame,
                        int32 compallsen);
int ptm_mgau_mllr_transform(ps_mgau_t *s,
                            ps_mllr_t *mllr);


#endif /*  __PTM_MGAU_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#if defined(__ADSPBLACKFIN__)
#elif !defined(_WIN32_WCE)
#include <sys/types.h>
#endif

/* SphinxBase headers */
#include <sphinx_config.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/fixpoint.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/bio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/prim_type.h>

/* Local headers */
#include "s2_semi_mgau.h"
#include "tied_mgau_common.h"

static ps_mgaufuncs_t s2_semi_mgau_funcs = {
    "s2_semi",
    s2_semi_mgau_frame_eval,      /* frame_eval */
    s2_semi_mgau_mllr_transform,  /* transform */
    s2_semi_mgau_free             /* free */
};

struct vqFeature_s {
    int32 score; /* score or distance */
    int32 codeword; /* codeword (vector index) */
};

static void
eval_topn_two(s2_semi_mgau_t *s, int32 feat, mfcc_t *z) //There are 2 eval topn functions defined.First one is this andisrenamed eval_topn_two
{
    int i, ceplen;
    vqFeature_t *topn;

    topn = s->f[feat];
    ceplen = s->veclen[feat];

    for (i = 0; i < s->max_topn; i++) {
        mfcc_t *mean, diff, sqdiff, compl; /* diff, diff^2, component likelihood */
        vqFeature_t vtmp;
        mfcc_t *var, d;
        mfcc_t *obs;
        int32 cw, j;

        cw = topn[i].codeword;
        mean = s->means[feat][0] + cw * ceplen;
        var = s->vars[feat][0] + cw * ceplen;
        d = s->dets[feat][cw];
        obs = z;
        for (j = 0; j < ceplen; j++) {
            diff = *obs++ - *mean++;
            sqdiff = MFCCMUL(diff, diff);
            compl = MFCCMUL(sqdiff, *var);
            d = GMMSUB(d, compl);
            ++var;
        }
        topn[i].score = (int32)d;
        if (i == 0)
            continue;
        vtmp = topn[i];
        for (j = i - 1; j >= 0 && (int32)d > topn[j].score; j--) {
            topn[j + 1] = topn[j];
        }
        topn[j + 1] = vtmp;
    }
}

static void
eval_cb_two(s2_semi_mgau_t *s, int32 feat, mfcc_t *z) //###########################################################################
{
    vqFeature_t *worst, *best, *topn;
    mfcc_t *mean;
    mfcc_t *var, *det, *detP, *detE;
    int32 i, ceplen;

    best = topn = s->f[feat];
    worst = topn + (s->max_topn - 1);
    mean = s->means[feat][0];
    var = s->vars[feat][0];
    det = s->dets[feat];
    detE = det + s->n_density;
    ceplen = s->veclen[feat];

    for (detP = det; detP < detE; ++detP) {
        mfcc_t diff, sqdiff, compl; /* diff, diff^2, component likelihood */
        mfcc_t d;
        mfcc_t *obs;
        vqFeature_t *cur;
        int32 cw, j;

        d = *detP;
        obs = z;
        cw = detP - det;
        for (j = 0; (j < ceplen) && (d >= worst->score); ++j) {
            diff = *obs++ - *mean++;
            sqdiff = MFCCMUL(diff, diff);
            compl = MFCCMUL(sqdiff, *var);
            d = GMMSUB(d, compl);
            ++var;
        }
        if (j < ceplen) {
            /* terminated early, so not in topn */
            mean += (ceplen - j);
            var += (ceplen - j);
            continue;
        }
        if ((int32)d < worst->score)
            continue;
        for (i = 0; i < s->max_topn; i++) {
            /* already there, so don't need to insert */
            if (topn[i].codeword == cw)
                break;
        }
        if (i < s->max_topn)
            continue;       /* already there.  Don't insert */
        /* remaining code inserts codeword and dist in correct spot */
        for (cur = worst - 1; cur >= best && (int32)d >= cur->score; --cur)
            memcpy(cur + 1, cur, sizeof(vqFeature_t));
        ++cur;
        cur->codeword = cw;
        cur->score = (int32)d;
    }
}

static void
mgau_dist(s2_semi_mgau_t * s, int32 frame, int32 feat, mfcc_t * z)
{
    eval_topn_two(s, feat, z); 

    /* If this frame is skipped, do nothing else. */
    if (frame % s->ds_ratio)
        return;

    /* Evaluate the rest of the codebook (or subset thereof). */
    eval_cb_two(s, feat, z);
}

static int
mgau_norm(s2_semi_mgau_t *s, int feat)
{
    int32 norm;
    int j;

    /* Compute quantized normalizing constant. */
    norm = s->f[feat][0].score >> SENSCR_SHIFT;

    /* Normalize the scores, negate them, and clamp their dynamic range. */
    for (j = 0; j < s->max_topn; ++j) {
        s->f[feat][j].score = -((s->f[feat][j].score >> SENSCR_SHIFT) - norm);
        if (s->f[feat][j].score > MAX_NEG_ASCR)
            s->f[feat][j].score = MAX_NEG_ASCR;
        if (s->topn_beam[feat] && s->f[feat][j].score > s->topn_beam[feat])
            break;
    }
    return j;
}

static int32
get_scores_8b_feat_6(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4, *pid_cw5;

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->mixw[i][s->f[i][4].codeword];
    pid_cw5 = s->mixw[i][s->f[i][5].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw3[sen] + s->f[i][3].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw4[sen] + s->f[i][4].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw5[sen] + s->f[i][5].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_5(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4;

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->mixw[i][s->f[i][4].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw3[sen] + s->f[i][3].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw4[sen] + s->f[i][4].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_4(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3;

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->mixw[i][s->f[i][3].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw3[sen] + s->f[i][3].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_3(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2;

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_2(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1;

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_1(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0;

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;
        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_any(s2_semi_mgau_t * s, int i, int topn,
                       int16 *senone_scores, uint8 *senone_active,
                       int32 n_senone_active)
{
    int32 j, k, l;

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        uint8 *pid_cw;
        int32 tmp;
        pid_cw = s->mixw[i][s->f[i][0].codeword];
        tmp = pid_cw[sen] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            pid_cw = s->mixw[i][s->f[i][k].codeword];
            tmp = fast_logmath_add(s->lmath_8b, tmp,
                                   pid_cw[sen] + s->f[i][k].score);
        }
        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat(s2_semi_mgau_t * s, int i, int topn,
                   int16 *senone_scores, uint8 *senone_active, int32 n_senone_active)
{
    switch (topn) {
    case 6:
        return get_scores_8b_feat_6(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 5:
        return get_scores_8b_feat_5(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 4:
        return get_scores_8b_feat_4(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 3:
        return get_scores_8b_feat_3(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 2:
        return get_scores_8b_feat_2(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 1:
        return get_scores_8b_feat_1(s, i, senone_scores,
                                    senone_active, n_senone_active);
    default:
        return get_scores_8b_feat_any(s, i, topn, senone_scores,
                                      senone_active, n_senone_active);
    }
}

static int32
get_scores_8b_feat_all(s2_semi_mgau_t * s, int i, int topn, int16 *senone_scores)
{
    int32 j, k;

    for (j = 0; j < s->n_sen; j++) {
        uint8 *pid_cw;
        int32 tmp;
        pid_cw = s->mixw[i][s->f[i][0].codeword];
        tmp = pid_cw[j] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            pid_cw = s->mixw[i][s->f[i][k].codeword];
            tmp = fast_logmath_add(s->lmath_8b, tmp,
                                   pid_cw[j] + s->f[i][k].score);
        }
        senone_scores[j] += tmp;
    }
    return 0;
}

static int32
get_scores_4b_feat_6(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4, *pid_cw5;
    uint8 w_den[6][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->mixw_cb[j] + s->f[i][2].score;
        w_den[3][j] = s->mixw_cb[j] + s->f[i][3].score;
        w_den[4][j] = s->mixw_cb[j] + s->f[i][4].score;
        w_den[5][j] = s->mixw_cb[j] + s->f[i][5].score;
    }

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->mixw[i][s->f[i][4].codeword];
    pid_cw5 = s->mixw[i][s->f[i][5].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
            cw = pid_cw5[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[5][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
            cw = pid_cw5[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[5][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_5(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4;
    uint8 w_den[5][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->mixw_cb[j] + s->f[i][2].score;
        w_den[3][j] = s->mixw_cb[j] + s->f[i][3].score;
        w_den[4][j] = s->mixw_cb[j] + s->f[i][4].score;
    }

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->mixw[i][s->f[i][4].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_4(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3;
    uint8 w_den[4][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->mixw_cb[j] + s->f[i][2].score;
        w_den[3][j] = s->mixw_cb[j] + s->f[i][3].score;
    }

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->mixw[i][s->f[i][3].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_3(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2;
    uint8 w_den[3][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->mixw_cb[j] + s->f[i][2].score;
    }

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->mixw[i][s->f[i][2].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_2(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1;
    uint8 w_den[2][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->mixw_cb[j] + s->f[i][1].score;
    }

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->mixw[i][s->f[i][1].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_1(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0;
    uint8 w_den[16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[j] = s->mixw_cb[j] + s->f[i][0].score;
    }

    pid_cw0 = s->mixw[i][s->f[i][0].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[cw];
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[cw];
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_any(s2_semi_mgau_t * s, int i, int topn,
                       int16 *senone_scores, uint8 *senone_active,
                       int32 n_senone_active)
{
    int32 j, k, l;

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;
        uint8 *pid_cw;
    
        pid_cw = s->mixw[i][s->f[i][0].codeword];
        if (n & 1)
            cw = pid_cw[n/2] >> 4;
        else
            cw = pid_cw[n/2] & 0x0f;
        tmp = s->mixw_cb[cw] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            pid_cw = s->mixw[i][s->f[i][k].codeword];
            if (n & 1)
                cw = pid_cw[n/2] >> 4;
            else
                cw = pid_cw[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp,
                                   s->mixw_cb[cw] + s->f[i][k].score);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat(s2_semi_mgau_t * s, int i, int topn,
                   int16 *senone_scores, uint8 *senone_active, int32 n_senone_active)
{
    switch (topn) {
    case 6:
        return get_scores_4b_feat_6(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 5:
        return get_scores_4b_feat_5(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 4:
        return get_scores_4b_feat_4(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 3:
        return get_scores_4b_feat_3(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 2:
        return get_scores_4b_feat_2(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 1:
        return get_scores_4b_feat_1(s, i, senone_scores,
                                    senone_active, n_senone_active);
    default:
        return get_scores_4b_feat_any(s, i, topn, senone_scores,
                                      senone_active, n_senone_active);
    }
}

static int32
get_scores_4b_feat_all(s2_semi_mgau_t * s, int i, int topn, int16 *senone_scores)
{
    int j, last_sen;

    j = 0;
    /* Number of senones is always even, but don't overrun if it isn't. */
    last_sen = s->n_sen & ~1;
    while (j < last_sen) {
        uint8 *pid_cw;
        int32 tmp0, tmp1;
        int k;

        pid_cw = s->mixw[i][s->f[i][0].codeword];
        tmp0 = s->mixw_cb[pid_cw[j/2] & 0x0f] + s->f[i][0].score;
        tmp1 = s->mixw_cb[pid_cw[j/2] >> 4] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            int32 w_den0, w_den1;

            pid_cw = s->mixw[i][s->f[i][k].codeword];
            w_den0 = s->mixw_cb[pid_cw[j/2] & 0x0f] + s->f[i][k].score;
            w_den1 = s->mixw_cb[pid_cw[j/2] >> 4] + s->f[i][k].score;
            tmp0 = fast_logmath_add(s->lmath_8b, tmp0, w_den0);
            tmp1 = fast_logmath_add(s->lmath_8b, tmp1, w_den1);
        }
        senone_scores[j++] += tmp0;
        senone_scores[j++] += tmp1;
    }
    return 0;
}

/*
 * Compute senone scores for the active senones.
 */
int32
s2_semi_mgau_frame_eval(ps_mgau_t *ps,
                        int16 *senone_scores,
                        uint8 *senone_active,
                        int32 n_senone_active,
			mfcc_t ** featbuf, int32 frame,
			int32 compallsen)
{
    s2_semi_mgau_t *s = (s2_semi_mgau_t *)ps;
    int i, topn_idx;

    memset(senone_scores, 0, s->n_sen * sizeof(*senone_scores));
    /* No bounds checking is done here, which just means you'll get
     * semi-random crap if you request a frame in the future or one
     * that's too far in the past. */
    topn_idx = frame % s->n_topn_hist;
    s->f = s->topn_hist[topn_idx];
    for (i = 0; i < s->n_feat; ++i) {
        /* For past frames this will already be computed. */
        if (frame >= ps_mgau_base(ps)->frame_idx) {
            vqFeature_t **lastf;
            if (topn_idx == 0)
                lastf = s->topn_hist[s->n_topn_hist-1];
            else
                lastf = s->topn_hist[topn_idx-1];
            memcpy(s->f[i], lastf[i], sizeof(vqFeature_t) * s->max_topn);
            mgau_dist(s, frame, i, featbuf[i]);
            s->topn_hist_n[topn_idx][i] = mgau_norm(s, i);
        }
        if (s->mixw_cb) {
            if (compallsen)
                get_scores_4b_feat_all(s, i, s->topn_hist_n[topn_idx][i], senone_scores);
            else
                get_scores_4b_feat(s, i, s->topn_hist_n[topn_idx][i], senone_scores,
                                   senone_active, n_senone_active);
        }
        else {
            if (compallsen)
                get_scores_8b_feat_all(s, i, s->topn_hist_n[topn_idx][i], senone_scores);
            else
                get_scores_8b_feat(s, i, s->topn_hist_n[topn_idx][i], senone_scores,
                                   senone_active, n_senone_active);
        }
    }

    return 0;
}

static int32
read_sendump_two(s2_semi_mgau_t *s, bin_mdef_t *mdef, char const *file)
{
    FILE *fp;
    char line[1000];
    int32 i, n, r, c;
    int32 do_swap, do_mmap;
    size_t offset;
    int n_clust = 0;
    int n_feat = s->n_feat;
    int n_density = s->n_density;
    int n_sen = bin_mdef_n_sen(mdef);
    int n_bits = 8;

    s->n_sen = n_sen; /* FIXME: Should have been done earlier*/ 
    do_mmap = cmd_ln_boolean_r(s->config, "-mmap");

    if ((fp = fopen(file, "rb")) == NULL)
        return -1;

    E_INFO("Loading senones from dump file %s\n", file);
    /* Read title size, title*/ 
    if (fread(&n, sizeof(int32), 1, fp) != 1) {
        E_ERROR_SYSTEM("Failed to read title size from %s", file);
        goto error_out;
    }
    /* This is extremely bogus */
    do_swap = 0;
    if (n < 1 || n > 999) {
        SWAP_INT32(&n);
        if (n < 1 || n > 999) {
            E_ERROR("Title length %x in dump file %s out of range\n", n, file);
            goto error_out;
        }
        do_swap = 1;
    }
    if (fread(line, sizeof(char), n, fp) != n) {
        E_ERROR_SYSTEM("Cannot read title");
        goto error_out;
    }
    if (line[n - 1] != '\0') {
        E_ERROR("Bad title in dump file\n");
        goto error_out;
    }
    E_INFO("%s\n", line);

    /* Read header size, header */
    if (fread(&n, sizeof(n), 1, fp) != 1) {
        E_ERROR_SYSTEM("Failed to read header size from %s", file);
        goto error_out;
    }
    if (do_swap) SWAP_INT32(&n);
    if (fread(line, sizeof(char), n, fp) != n) {
        E_ERROR_SYSTEM("Cannot read header");
        goto error_out;
    }
    if (line[n - 1] != '\0') {
        E_ERROR("Bad header in dump file\n");
        goto error_out;
    }

    /* Read other header strings until string length = 0*/ 
    for (;;) {
        if (fread(&n, sizeof(n), 1, fp) != 1) {
            E_ERROR_SYSTEM("Failed to read header string size from %s", file);
            goto error_out;
        }
        if (do_swap) SWAP_INT32(&n);
        if (n == 0)
            break;
        if (fread(line, sizeof(char), n, fp) != n) {
            E_ERROR_SYSTEM("Cannot read header");
            goto error_out;
        }
        /* Look for a cluster count, if present */
        if (!strncmp(line, "feature_count ", strlen("feature_count "))) {
            n_feat = atoi(line + strlen("feature_count "));
        }
        if (!strncmp(line, "mixture_count ", strlen("mixture_count "))) {
            n_density = atoi(line + strlen("mixture_count "));
        }
        if (!strncmp(line, "model_count ", strlen("model_count "))) {
            n_sen = atoi(line + strlen("model_count "));
        }
        if (!strncmp(line, "cluster_count ", strlen("cluster_count "))) {
            n_clust = atoi(line + strlen("cluster_count "));
        }
        if (!strncmp(line, "cluster_bits ", strlen("cluster_bits "))) {
            n_bits = atoi(line + strlen("cluster_bits "));
        }
    }

    /* Defaults for #rows, #columns in mixw array.*/ 
    c = n_sen;
    r = n_density;
    if (n_clust == 0) {
        /* Older mixw files have them here, and they might be padded.*/ 
        if (fread(&r, sizeof(r), 1, fp) != 1) {
            E_ERROR_SYSTEM("Cannot read #rows");
            goto error_out;
        }
        if (do_swap) SWAP_INT32(&r);
        if (fread(&c, sizeof(c), 1, fp) != 1) {
            E_ERROR_SYSTEM("Cannot read #columns");
            goto error_out;
        }
        if (do_swap) SWAP_INT32(&c);
        E_INFO("Rows: %d, Columns: %d\n", r, c);
    }

    if (n_feat != s->n_feat) {
        E_ERROR("Number of feature streams mismatch: %d != %d\n",
                n_feat, s->n_feat);
        goto error_out;
    }
    if (n_density != s->n_density) {
        E_ERROR("Number of densities mismatch: %d != %d\n",
                n_density, s->n_density);
        goto error_out;
    }
    if (n_sen != s->n_sen) {
        E_ERROR("Number of senones mismatch: %d != %d\n",
                n_sen, s->n_sen);
        goto error_out;
    }

    if (!((n_clust == 0) || (n_clust == 15) || (n_clust == 16))) {
        E_ERROR("Cluster count must be 0, 15, or 16\n");
        goto error_out;
    }
    if (n_clust == 15)
        ++n_clust;

    if (!((n_bits == 8) || (n_bits == 4))) {
        E_ERROR("Cluster count must be 4 or 8\n");
        goto error_out;
    }

    if (do_mmap) {
            E_INFO("Using memory-mapped I/O for senones\n");
    }
    offset = ftell(fp);

    /* Allocate memory for pdfs (or memory map them) */
    if (do_mmap) {
        s->sendump_mmap = mmio_file_read(file);
        /* Get cluster codebook if any. */
        if (n_clust) {
            s->mixw_cb = ((uint8 *) mmio_file_ptr(s->sendump_mmap)) + offset;
            offset += n_clust;
        }
    }
    else {
        /* Get cluster codebook if any. */
        if (n_clust) {
            s->mixw_cb = ckd_calloc(1, n_clust);
            if (fread(s->mixw_cb, 1, n_clust, fp) != (size_t) n_clust) {
                E_ERROR("Failed to read %d bytes from sendump\n", n_clust);
                goto error_out;
            }
        }
    }

    /* Set up pointers, or read, or whatever*/ 
    if (s->sendump_mmap) {
        s->mixw = ckd_calloc_2d(s->n_feat, n_density, sizeof(*s->mixw));
        for (n = 0; n < n_feat; n++) {
            int step = c;
            if (n_bits == 4)
                step = (step + 1) / 2;
            for (i = 0; i < r; i++) {
                s->mixw[n][i] = ((uint8 *) mmio_file_ptr(s->sendump_mmap)) + offset;
                offset += step;
            }
        }
    }
    else {
        s->mixw = ckd_calloc_3d(n_feat, n_density, n_sen, sizeof(***s->mixw));
        /* Read pdf values and ids */
        for (n = 0; n < n_feat; n++) {
            int step = c;
            if (n_bits == 4)
                step = (step + 1) / 2;
            for (i = 0; i < r; i++) {
                if (fread(s->mixw[n][i], sizeof(***s->mixw), step, fp)
                    != (size_t) step) {
                    E_ERROR("Failed to read %d bytes from sendump\n", step);
                    goto error_out;
                }
            }
        }
    }

    fclose(fp);
    return 0;
error_out:
    fclose(fp);
    return -1;
}

static int32
read_mixw_two(s2_semi_mgau_t * s, char const *file_name, double SmoothMin)
{
    char **argname, **argval;
    char eofchk;
    FILE *fp;
    int32 byteswap, chksum_present;
    uint32 chksum;
    float32 *pdf;
    int32 i, f, c, n;
    int32 n_sen;
    int32 n_feat;
    int32 n_comp;
    int32 n_err;

    E_INFO("Reading mixture weights file '%s'\n", file_name);

    if ((fp = fopen(file_name, "rb")) == NULL)
        E_FATAL_SYSTEM("Failed to open mixture weights file '%s' for reading", file_name);

    /* Read header, including argument-value info and 32-bit byteorder magic*/ 
    if (bio_readhdr(fp, &argname, &argval, &byteswap) < 0)
        E_FATAL("Failed to read header from file '%s'\n", file_name);

    /* Parse argument-value list */
    chksum_present = 0;
    for (i = 0; argname[i]; i++) {
        if (strcmp(argname[i], "version") == 0) {
            if (strcmp(argval[i], MGAU_MIXW_VERSION) != 0)
                E_WARN("Version mismatch(%s): %s, expecting %s\n",
                       file_name, argval[i], MGAU_MIXW_VERSION);
        }
        else if (strcmp(argname[i], "chksum0") == 0) {
            chksum_present = 1; /* Ignore the associated value*/ 
        }
    }
    bio_hdrarg_free(argname, argval);
    argname = argval = NULL;

    chksum = 0;

    /* Read #senones, #features, #codewords, arraysize */
    if ((bio_fread(&n_sen, sizeof(int32), 1, fp, byteswap, &chksum) != 1)
        || (bio_fread(&n_feat, sizeof(int32), 1, fp, byteswap, &chksum) !=
            1)
        || (bio_fread(&n_comp, sizeof(int32), 1, fp, byteswap, &chksum) !=
            1)
        || (bio_fread(&n, sizeof(int32), 1, fp, byteswap, &chksum) != 1)) {
        E_FATAL("bio_fread(%s) (arraysize) failed\n", file_name);
    }
    if (n_feat != s->n_feat)
        E_FATAL("#Features streams(%d) != %d\n", n_feat, s->n_feat);
    if (n != n_sen * n_feat * n_comp) {
        E_FATAL
            ("%s: #float32s(%d) doesn't match header dimensions: %d x %d x %d\n",
             file_name, i, n_sen, n_feat, n_comp);
    }

    /* n_sen = number of mixture weights per codeword, which is
     * fixed at the number of senones since we have only one codebook.
     */
    s->n_sen = n_sen;

    /* Quantized mixture weight arrays.*/ 
    s->mixw = ckd_calloc_3d(s->n_feat, s->n_density, n_sen, sizeof(***s->mixw));

    /* Temporary structure to read in floats before conversion to (int32) logs3*/ 
    pdf = (float32 *) ckd_calloc(n_comp, sizeof(float32));

    /* Read senone probs data, normalize, floor, convert to logs3, truncate to 8 bits*/ 
    n_err = 0;
    for (i = 0; i < n_sen; i++) {
        for (f = 0; f < n_feat; f++) {
            if (bio_fread((void *) pdf, sizeof(float32),
                          n_comp, fp, byteswap, &chksum) != n_comp) {
                E_FATAL("bio_fread(%s) (arraydata) failed\n", file_name);
            }

            /* Normalize and floor */
            if (vector_sum_norm(pdf, n_comp) <= 0.0)
                n_err++;
            vector_floor(pdf, n_comp, SmoothMin);
            vector_sum_norm(pdf, n_comp);

            /* Convert to LOG, quantize, and transpose*/ 
            for (c = 0; c < n_comp; c++) {
                int32 qscr;

                qscr = -logmath_log(s->lmath_8b, pdf[c]);
                if ((qscr > MAX_NEG_MIXW) || (qscr < 0))
                    qscr = MAX_NEG_MIXW;
                s->mixw[f][c][i] = qscr;
            }
        }
    }
    if (n_err > 0)
        E_WARN("Weight normalization failed for %d mixture weights components\n", n_err);

    ckd_free(pdf);

    if (chksum_present)
        bio_verify_chksum(fp, byteswap, chksum);

    if (fread(&eofchk, 1, 1, fp) == 1)
        E_FATAL("More data than expected in %s\n", file_name);

    fclose(fp);

    E_INFO("Read %d x %d x %d mixture weights\n", n_sen, n_feat, n_comp);
    return n_sen;
}

static int
split_topn(char const *str, uint8 *out, int nfeat)
{
    char *topn_list = ckd_salloc(str);
    char *c, *cc;
    int i, maxn;

    c = topn_list;
    i = 0;
    maxn = 0;
    while (i < nfeat && (cc = strchr(c, ',')) != NULL) {
        *cc = '\0';
        out[i] = atoi(c);
        if (out[i] > maxn) maxn = out[i];
        c = cc + 1;
        ++i;
    }
    if (i < nfeat && *c != '\0') {
        out[i] = atoi(c);
        if (out[i] > maxn) maxn = out[i];
        ++i;
    }
    while (i < nfeat)
        out[i++] = maxn;

    ckd_free(topn_list);
    return maxn;
}


ps_mgau_t *
s2_semi_mgau_init(acmod_t *acmod)
{
    s2_semi_mgau_t *s;
    ps_mgau_t *ps;
    char const *sendump_path;
    int i;

    s = ckd_calloc(1, sizeof(*s));
    s->config = acmod->config;

    s->lmath = logmath_retain(acmod->lmath);
    /* Log-add table. */
    s->lmath_8b = logmath_init(logmath_get_base(acmod->lmath), SENSCR_SHIFT, TRUE);
    if (s->lmath_8b == NULL)
        goto error_out;
    /* Ensure that it is only 8 bits wide so that fast_logmath_add() works. */
    if (logmath_get_width(s->lmath_8b) != 1) {
        E_ERROR("Log base %f is too small to represent add table in 8 bits\n",
                logmath_get_base(s->lmath_8b));
        goto error_out;
    }

    /* Read means and variances. */
    if ((s->g = gauden_init(cmd_ln_str_r(s->config, "-mean"),
                            cmd_ln_str_r(s->config, "-var"),
                            cmd_ln_float32_r(s->config, "-varfloor"),
                            s->lmath)) == NULL)
        goto error_out;
    /* Currently only a single codebook is supported. */
    if (s->g->n_mgau != 1)
        goto error_out;
    /* FIXME: maintaining pointers for convenience for now */
    s->means = s->g->mean[0];
    s->vars = s->g->var[0];
    s->dets = s->g->det[0];
    s->veclen = s->g->featlen;    
    /* Verify n_feat and veclen, against acmod. */
    s->n_feat = s->g->n_feat;
    if (s->n_feat != feat_dimension1(acmod->fcb)) {
        E_ERROR("Number of streams does not match: %d != %d\n",
                s->n_feat, feat_dimension1(acmod->fcb));
        goto error_out;
    }
    for (i = 0; i < s->n_feat; ++i) {
        if (s->veclen[i] != feat_dimension2(acmod->fcb, i)) {
            E_ERROR("Dimension of stream %d does not match: %d != %d\n",
                    i, s->veclen[i], feat_dimension2(acmod->fcb, i));
            goto error_out;
        }
    }
    s->n_density = s->g->n_density;
    /* Read mixture weights */
    if ((sendump_path = cmd_ln_str_r(s->config, "-sendump"))) {
        if (read_sendump_two(s, acmod->mdef, sendump_path) < 0) {
            goto error_out;
        }
    }
    else {
        if (read_mixw_two(s, cmd_ln_str_r(s->config, "-mixw"),
                      cmd_ln_float32_r(s->config, "-mixwfloor")) < 0) {
            goto error_out;
        }
    }
    s->ds_ratio = cmd_ln_int32_r(s->config, "-ds");

    /* Determine top-N for each feature */
    s->topn_beam = ckd_calloc(s->n_feat, sizeof(*s->topn_beam));
    s->max_topn = cmd_ln_int32_r(s->config, "-topn");
    split_topn(cmd_ln_str_r(s->config, "-topn_beam"), s->topn_beam, s->n_feat);
    E_INFO("Maximum top-N: %d ", s->max_topn);
    E_INFOCONT("Top-N beams:");
    for (i = 0; i < s->n_feat; ++i) {
        E_INFOCONT(" %d", s->topn_beam[i]);
    }
    E_INFOCONT("\n");

    /* Top-N scores from recent frames */
    s->n_topn_hist = cmd_ln_int32_r(s->config, "-pl_window") + 2;
    s->topn_hist = (vqFeature_t ***)
        ckd_calloc_3d(s->n_topn_hist, s->n_feat, s->max_topn,
                      sizeof(***s->topn_hist));
    s->topn_hist_n = ckd_calloc_2d(s->n_topn_hist, s->n_feat,
                                   sizeof(**s->topn_hist_n));
    for (i = 0; i < s->n_topn_hist; ++i) {
        int j;
        for (j = 0; j < s->n_feat; ++j) {
            int k;
            for (k = 0; k < s->max_topn; ++k) {
                s->topn_hist[i][j][k].score = WORST_DIST;
                s->topn_hist[i][j][k].codeword = k;
            }
        }
    }

    ps = (ps_mgau_t *)s;
    ps->vt = &s2_semi_mgau_funcs;
    return ps;
error_out:
    s2_semi_mgau_free(ps_mgau_base(s));
    return NULL;
}

int
s2_semi_mgau_mllr_transform(ps_mgau_t *ps,
                            ps_mllr_t *mllr)
{
    s2_semi_mgau_t *s = (s2_semi_mgau_t *)ps;
    return gauden_mllr_transform(s->g, mllr, s->config);
}

void
s2_semi_mgau_free(ps_mgau_t *ps)
{
    s2_semi_mgau_t *s = (s2_semi_mgau_t *)ps;

    logmath_free(s->lmath);
    logmath_free(s->lmath_8b);
    if (s->sendump_mmap) {
        ckd_free_2d(s->mixw); 
        mmio_file_unmap(s->sendump_mmap);
    }
    else {
        ckd_free_3d(s->mixw);
        if (s->mixw_cb)
            ckd_free(s->mixw_cb);
    }
    gauden_free(s->g);
    ckd_free(s->topn_beam);
    ckd_free_2d(s->topn_hist_n);
    ckd_free_3d((void **)s->topn_hist);
    ckd_free(s);
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * Interface for "semi-continuous vector quantization", a.k.a. Sphinx2
 * fast GMM computation.
 */

#ifndef __S2_SEMI_MGAU_H__
#define __S2_SEMI_MGAU_H__

/* SphinxBase headesr. */
#include <sphinxbase/fe.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/mmio.h>

/* Local headers. */
#include "acmod.h"
#include "hmm.h"
#include "bin_mdef.h"
#include "ms_gauden.h"

typedef struct vqFeature_s vqFeature_t;

typedef struct s2_semi_mgau_s s2_semi_mgau_t;
struct s2_semi_mgau_s {
    ps_mgau_t base;     /**< base structure. */
    cmd_ln_t *config;   /* configuration parameters */

    gauden_t *g;        /* Set of Gaussians (pointers below point in here and will go away soon) */
    mfcc_t  ***means;	/* mean vectors foreach feature, density */
    mfcc_t  ***vars;	/* inverse var vectors foreach feature, density */
    mfcc_t  **dets;	/* det values foreach cb, feature */

    uint8 ***mixw;     /* mixture weight distributions */
    mmio_file_t *sendump_mmap;/* memory map for mixw (or NULL if not mmap) */

    uint8 *mixw_cb;    /* mixture weight codebook, if any (assume it contains 16 values) */
    int32 *veclen;	/* Length of feature streams */
    int16 n_feat;	/* Number of feature streams */
    int16 n_density;	/* Number of mixtures per codebook */
    int32 n_sen;	/* Number of senones */
    uint8 *topn_beam;   /* Beam for determining per-frame top-N densities */
    int16 max_topn;
    int16 ds_ratio;

    vqFeature_t ***topn_hist; /**< Top-N scores and codewords for past frames. */
    uint8 **topn_hist_n;      /**< Variable top-N for past frames. */
    vqFeature_t **f;          /**< Topn-N for currently scoring frame. */
    int n_topn_hist;          /**< Number of past frames tracked. */

    /* Log-add table for compressed values. */
    logmath_t *lmath_8b;
    /* Log-add object for reloading means/variances. */
    logmath_t *lmath;
};

ps_mgau_t *s2_semi_mgau_init(acmod_t *acmod);
void s2_semi_mgau_free(ps_mgau_t *s);
int s2_semi_mgau_frame_eval(ps_mgau_t *s,
                            int16 *senone_scores,
                            uint8 *senone_active,
                            int32 n_senone_active,
                            mfcc_t **featbuf,
                            int32 frame,
                            int32 compallsen);
int s2_semi_mgau_mllr_transform(ps_mgau_t *s,
                                ps_mllr_t *mllr);


#endif /*  __S2_SEMI_MGAU_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * s3types.h -- Types specific to s3 decoder.
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1999 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * $Log: s3types.h,v $
 * Revision 1.16  2006/02/22 19:57:57  arthchan2003
 * Merged from branch SPHINX3_5_2_RCI_IRII_BRANCH: Increase the size of MAX_S3CIPID from 127 to 32767.  This will make Chinese Mandarin setup works.
 *
 * Revision 1.15.4.1  2005/10/09 19:53:09  arthchan2003
 * Changed the maximum number of CI PID from 127 to 32767, this will allow us to take care of Chinese syllable, Chinese initial/final and even Cantononese.  It might still cause us problem in Turkish.
 *
 * Revision 1.15  2005/06/21 20:54:44  arthchan2003
 * 1, Added $ keyword. 2, make a small change for compilation purpose.
 *
 * Revision 1.5  2005/06/16 04:59:09  archan
 * Sphinx3 to s3.generic, a gentle-refactored version of Dave's change in senone scale.
 *
 * Revision 1.4  2005/06/15 21:48:56  archan
 * Sphinx3 to s3.generic: Changed noinst_HEADERS to pkginclude_HEADERS.  This make all the headers to be installed.
 *
 * Revision 1.3  2005/03/30 01:22:47  archan
 * Fixed mistakes in last updates. Add
 *
 * 
 * 13-May-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Changed typedef source for s3ssid_t from int32 to s3pid_t.
 * 		Changed s3senid_t from int16 to int32 (to conform with composite senid
 * 		which is int32).
 * 
 * 04-May-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added senone sequence ID (s3ssid_t).
 * 
 * 12-Jul-95	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Started.
 */


#ifndef _S3_S3TYPES_H_
#define _S3_S3TYPES_H_

#include <float.h>
#include <assert.h>

#include <sphinxbase/prim_type.h>
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>

/** \file s3types.h
 * \brief Size definition of semantically units. Common for both s3 and s3.X decoder. 
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

/**
 * Size definitions for more semantially meaningful units.
 * Illegal value definitions, limits, and tests for specific types.
 * NOTE: Types will be either int32 or smaller; only smaller ones may be unsigned (i.e.,
 * no type will be uint32).
 */

typedef int16		s3cipid_t;	/** Ci phone id */
#define BAD_S3CIPID	((s3cipid_t) -1)
#define NOT_S3CIPID(p)	((p)<0)
#define IS_S3CIPID(p)	((p)>=0)
#define MAX_S3CIPID	32767

/*#define MAX_S3CIPID	127*/

typedef int32		s3pid_t;	/** Phone id (triphone or ciphone) */
#define BAD_S3PID	((s3pid_t) -1)
#define NOT_S3PID(p)	((p)<0)
#define IS_S3PID(p)	((p)>=0)
#define MAX_S3PID	((int32)0x7ffffffe)

typedef uint16		s3ssid_t;	/** Senone sequence id (triphone or ciphone) */
#define BAD_S3SSID	((s3ssid_t) 0xffff)
#define NOT_S3SSID(p)	((p) == BAD_S3SSID)
#define IS_S3SSID(p)	((p) != BAD_S3SSID)
#define MAX_S3SSID	((s3ssid_t)0xfffe)

typedef int32		s3tmatid_t;	/** Transition matrix id; there can be as many as pids */
#define BAD_S3TMATID	((s3tmatid_t) -1)
#define NOT_S3TMATID(t)	((t)<0)
#define IS_S3TMATID(t)	((t)>=0)
#define MAX_S3TMATID	((int32)0x7ffffffe)

typedef int32		s3wid_t;	/** Dictionary word id */
#define BAD_S3WID	((s3wid_t) -1)
#define NOT_S3WID(w)	((w)<0)
#define IS_S3WID(w)	((w)>=0)
#define MAX_S3WID	((int32)0x7ffffffe)

typedef uint16		s3lmwid_t;	/** LM word id (uint16 for conserving space) */
#define BAD_S3LMWID	((s3lmwid_t) 0xffff)
#define NOT_S3LMWID(w)	((w)==BAD_S3LMWID)
#define IS_S3LMWID(w)	((w)!=BAD_S3LMWID)
#define MAX_S3LMWID	((uint32)0xfffe)
#define BAD_LMCLASSID   (-1)

typedef uint32		s3lmwid32_t;	/** LM word id (uint32 for conserving space) */
#define BAD_S3LMWID32	((s3lmwid32_t) 0x0fffffff)
#define NOT_S3LMWID32(w)  ((w)==BAD_S3LMWID32)
#define IS_S3LMWID32(w)	((w)!=BAD_S3LMWID32)
#define MAX_S3LMWID32	((uint32)0xfffffffe)

/* Generic macro that is applicable to both uint16 and uint32
   Careful with efficiency issue. 

   Also, please don't use BAD_S3LATID(l);
*/

#define BAD_LMWID(lm)      (lm->is32bits? BAD_S3LMWID32 : BAD_S3LMWID)
#define NOT_LMWID(lm,w)    (lm->is32bits? NOT_S3LMWID32(w): NOT_S3LMWID(w))
#define IS_LMWID(lm,w)     (lm->is32bits? IS_S3LMWID32(w): IS_S3LMWID(w))
#define MAX_LMWID(lm)      (lm->is32bits? MAX_S3LMWID32: MAX_S3LMWID)

typedef int32		s3latid_t;	/** Lattice entry id */
#define BAD_S3LATID	((s3latid_t) -1)
#define NOT_S3LATID(l)	((l)<0)
#define IS_S3LATID(l)	((l)>=0)
#define MAX_S3LATID	((int32)0x7ffffffe)

typedef int16   	s3frmid_t;	/** Frame id (must be SIGNED integer) */
#define BAD_S3FRMID	((s3frmid_t) -1)
#define NOT_S3FRMID(f)	((f)<0)
#define IS_S3FRMID(f)	((f)>=0)
#define MAX_S3FRMID	((int32)0x7ffe)

typedef uint16   	s3senid_t;	/** Senone id */
#define BAD_S3SENID	((s3senid_t) 0xffff)
#define NOT_S3SENID(s)	((s) == BAD_S3SENID)
#define IS_S3SENID(s)	((s) != BAD_S3SENID)
#define MAX_S3SENID	((int16)0x7ffe)

typedef int16   	s3mgauid_t;	/** Mixture-gaussian codebook id */
#define BAD_S3MGAUID	((s3mgauid_t) -1)
#define NOT_S3MGAUID(m)	((m)<0)
#define IS_S3MGAUID(m)	((m)>=0)
#define MAX_S3MGAUID	((int32)0x00007ffe)


#define S3_LOGPROB_ZERO		((int32) 0xc8000000)	/** Integer version of log of zero Approx -infinity!! */
#define S3_LOGPROB_ZERO_F	((float32) -1e30)	/** Float version of log of zero Approx -infinity!! */

#define RENORM_THRESH     ((int32) ((S3_LOGPROB_ZERO)>>1))       /** Bestscore getting close to 0 */

#define S3_SUCCESS      0
#define S3_ERROR        -1
#define S3_WARNING      -2

/** The maximum # of states for any given acoustic model */
#define MAX_N_STATE     20

/** The maximum # of attributes associated with any
 * given acoustic model */
#define MAX_N_ATTRIB    5

#ifndef TRUE
#define TRUE  1
#define FALSE 0 /* assume that true is never defined w/o false */
#endif

/* Timer for elapsed I/O time */
#define IO_ELAPSED      0

/* Timer for utt processing elapsed time */
#define UTT_ELAPSED     1
#define UTT_IO_ELAPSED  2
#define UTT_BW_ELAPSED  3

#define TYING_NON_EMITTING      (0xffffffff)
#define TYING_NO_ID             (0xffffffff)

#define MAX_VERSION_LEN 128

#define MEG *1024*1024

#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif

#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file state_align_search.c State (and phone and word) alignment search.
 */

#include "state_align_search.h"

static int
state_align_search_start(ps_search_t *search)
{
    state_align_search_t *sas = (state_align_search_t *)search;

    /* Activate the initial state. */
    hmm_enter(sas->hmms, 0, 0, 0);

    return 0;
}

static void
renormalize_hmms_two(state_align_search_t *sas, int frame_idx, int32 norm)
{
    int i;
    for (i = 0; i < sas->n_phones; ++i)
        hmm_normalize(sas->hmms + i, norm);
}

static int32
evaluate_hmms_two(state_align_search_t *sas, int16 const *senscr, int frame_idx)
{
    int32 bs = WORST_SCORE;
    int i;

    hmm_context_set_senscore(sas->hmmctx, senscr);

    for (i = 0; i < sas->n_phones; ++i) {
        hmm_t *hmm = sas->hmms + i;
        int32 score;

        if (hmm_frame(hmm) < frame_idx)
            continue;
        score = hmm_vit_eval(hmm);
        if (score BETTER_THAN bs) {
            bs = score;
        }
    }
    return bs;
}

static void
prune_hmms_two(state_align_search_t *sas, int frame_idx)
{
    int nf = frame_idx + 1;
    int i;

    /* Check all phones to see if they remain active in the next frame. */
    for (i = 0; i < sas->n_phones; ++i) {
        hmm_t *hmm = sas->hmms + i;
        if (hmm_frame(hmm) < frame_idx)
            continue;
        hmm_frame(hmm) = nf;
    }
}

static void
phone_transition_two(state_align_search_t *sas, int frame_idx)
{
    int nf = frame_idx + 1;
    int i;

    for (i = 0; i < sas->n_phones - 1; ++i) {
        hmm_t *hmm, *nhmm;
        int32 newphone_score;

        hmm = sas->hmms + i;
        if (hmm_frame(hmm) != nf)
            continue;

        newphone_score = hmm_out_score(hmm);
        /* Transition into next phone using the usual Viterbi rule. */
        nhmm = hmm + 1;
        if (hmm_frame(nhmm) < frame_idx
            || newphone_score BETTER_THAN hmm_in_score(nhmm)) {
            hmm_enter(nhmm, newphone_score, hmm_out_history(hmm), nf);
        }
    }
}

#define TOKEN_STEP 20
static void
extend_tokenstack(state_align_search_t *sas, int frame_idx)
{
    if (frame_idx >= sas->n_fr_alloc) {
        sas->n_fr_alloc = frame_idx + TOKEN_STEP + 1;
        sas->tokens = ckd_realloc(sas->tokens,
                                  sas->n_emit_state * sas->n_fr_alloc
                                  * sizeof(*sas->tokens));
    }
    memset(sas->tokens + frame_idx * sas->n_emit_state, 0xff,
           sas->n_emit_state * sizeof(*sas->tokens));
}

static void
record_transitions(state_align_search_t *sas, int frame_idx)
{
    uint16 *tokens;
    int i;

    /* Push another frame of tokens on the stack. */
    extend_tokenstack(sas, frame_idx);
    tokens = sas->tokens + frame_idx * sas->n_emit_state;

    /* Scan all active HMMs */
    for (i = 0; i < sas->n_phones; ++i) {
        hmm_t *hmm = sas->hmms + i;
        int j;

        if (hmm_frame(hmm) < frame_idx)
            continue;
        for (j = 0; j < sas->hmmctx->n_emit_state; ++j) {
            int state_idx = i * sas->hmmctx->n_emit_state + j;
            /* Record their backpointers on the token stack. */
            tokens[state_idx] = hmm_history(hmm, j);
            /* Update backpointer fields with state index. */
            hmm_history(hmm, j) = state_idx;
        }
    }
}

static int
state_align_search_step(ps_search_t *search, int frame_idx)
{
    state_align_search_t *sas = (state_align_search_t *)search;
    acmod_t *acmod = ps_search_acmod(search);
    int16 const *senscr;
    int i;

    /* Calculate senone scores. */
    for (i = 0; i < sas->n_phones; ++i)
        acmod_activate_hmm(acmod, sas->hmms + i);
    senscr = acmod_score(acmod, &frame_idx);

    /* Renormalize here if needed. */
    /* FIXME: Make sure to (unit-)test this!!! */
    if ((sas->best_score - 0x300000) WORSE_THAN WORST_SCORE) {
        E_INFO("Renormalizing Scores at frame %d, best score %d\n",
               frame_idx, sas->best_score);
        renormalize_hmms_two(sas, frame_idx, sas->best_score);
    }
    
    /* Viterbi step. */
    sas->best_score = evaluate_hmms_two(sas, senscr, frame_idx);
    prune_hmms_two(sas, frame_idx);

    /* Transition out of non-emitting states. */
    phone_transition_two(sas, frame_idx);

    /* Generate new tokens from best path results. */
    record_transitions(sas, frame_idx);

    /* Update frame counter */
    sas->frame = frame_idx;

    return 0;
}

static int
state_align_search_finish(ps_search_t *search)
{
    state_align_search_t *sas = (state_align_search_t *)search;
    hmm_t *final_phone = sas->hmms + sas->n_phones - 1;
    ps_alignment_iter_t *itor;
    ps_alignment_entry_t *ent;
    int next_state, next_start, state, frame;

    /* Best state exiting the last frame. */
    next_state = state = hmm_out_history(final_phone);
    if (state == 0xffff) {
        E_ERROR("Failed to reach final state in alignment\n");
        return -1;
    }
    itor = ps_alignment_states(sas->al);
    next_start = sas->frame + 1;
    for (frame = sas->frame - 1; frame >= 0; --frame) {
        state = sas->tokens[frame * sas->n_emit_state + state];
        /* State boundary, update alignment entry for next state. */
        if (state != next_state) {
            itor = ps_alignment_iter_goto(itor, next_state);
            assert(itor != NULL);
            ent = ps_alignment_iter_get(itor);
            ent->start = frame + 1;
            ent->duration = next_start - ent->start;
            E_DEBUG(1,("state %d start %d end %d\n", next_state,
                       ent->start, next_start));
            next_state = state;
            next_start = frame + 1;
        }
    }
    /* Update alignment entry for initial state. */
    itor = ps_alignment_iter_goto(itor, 0);
    assert(itor != NULL);
    ent = ps_alignment_iter_get(itor);
    ent->start = 0;
    ent->duration = next_start;
    E_DEBUG(1,("state %d start %d end %d\n", 0,
               ent->start, next_start));
    ps_alignment_iter_free(itor);
    ps_alignment_propagate(sas->al);

    return 0;
}

static int
state_align_search_reinit(ps_search_t *search, dict_t *dict, dict2pid_t *d2p)
{
    /* This does nothing. */
    return 0;
}

static void
state_align_search_free(ps_search_t *search)
{
    state_align_search_t *sas = (state_align_search_t *)search;
    ps_search_deinit(search);
    ckd_free(sas->hmms);
    ckd_free(sas->tokens);
    hmm_context_free(sas->hmmctx);
    ckd_free(sas);
}

static ps_searchfuncs_t state_align_search_funcs = {
    /* name: */   "state_align",
    /* start: */  state_align_search_start,
    /* step: */   state_align_search_step,
    /* finish: */ state_align_search_finish,
    /* reinit: */ state_align_search_reinit,
    /* free: */   state_align_search_free,
    /* lattice: */  NULL,
    /* hyp: */      NULL,
    /* prob: */     NULL,
    /* seg_iter: */ NULL,
};

ps_search_t *
state_align_search_init(cmd_ln_t *config,
                        acmod_t *acmod,
                        ps_alignment_t *al)
{
    state_align_search_t *sas;
    ps_alignment_iter_t *itor;
    hmm_t *hmm;

    sas = ckd_calloc(1, sizeof(*sas));
    ps_search_init(ps_search_base(sas), &state_align_search_funcs,
                   config, acmod, al->d2p->dict, al->d2p);
    sas->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
                                   acmod->tmat->tp, NULL, acmod->mdef->sseq);
    if (sas->hmmctx == NULL) {
        ckd_free(sas);
        return NULL;
    }
    sas->al = al;

    /* Generate HMM vector from phone level of alignment. */
    sas->n_phones = ps_alignment_n_phones(al);
    sas->n_emit_state = ps_alignment_n_states(al);
    sas->hmms = ckd_calloc(sas->n_phones, sizeof(*sas->hmms));
    for (hmm = sas->hmms, itor = ps_alignment_phones(al); itor;
         ++hmm, itor = ps_alignment_iter_next(itor)) {
        ps_alignment_entry_t *ent = ps_alignment_iter_get(itor);
        hmm_init(sas->hmmctx, hmm, FALSE,
                 ent->id.pid.ssid, ent->id.pid.tmatid);
    }
    return ps_search_base(sas);
}
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file state_align_search.h State (and phone and word) alignment search.
 */

#ifndef __STATE_ALIGN_SEARCH_H__
#define __STATE_ALIGN_SEARCH_H__

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>

/* Local headers. */
#include "pocketsphinx_internal.h"
#include "ps_alignment.h"
#include "hmm.h"

/**
 * Phone loop search structure.
 */
struct state_align_search_s {
    ps_search_t base;       /**< Base search structure. */
    hmm_context_t *hmmctx;  /**< HMM context structure. */
    ps_alignment_t *al;     /**< Alignment structure being operated on. */
    hmm_t *hmms;            /**< Vector of HMMs corresponding to phone level. */
    int n_phones;	    /**< Number of HMMs (phones). */

    int frame;              /**< Current frame being processed. */
    int32 best_score;       /**< Best score in current frame. */

    int n_emit_state;       /**< Number of emitting states (tokens per frame) */
    uint16 *tokens;         /**< Tokens (backpointers) for state alignment. */
    int n_fr_alloc;         /**< Number of frames of tokens allocated. */
};
typedef struct state_align_search_s state_align_search_t;

ps_search_t *state_align_search_init(cmd_ln_t *config,
                                     acmod_t *acmod,
                                     ps_alignment_t *al);

#endif /* __STATE_ALIGN_SEARCH_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file tied_mgau_common.h
 * @brief Common code shared between SC and PTM (tied-state) models.
 */

#ifndef __TIED_MGAU_COMMON_H__
#define __TIED_MGAU_COMMON_H__

#include <sphinxbase/logmath.h>
#include <sphinxbase/fixpoint.h>

#define MGAU_MIXW_VERSION	"1.0"   /* Sphinx-3 file format version for mixw */
#define MGAU_PARAM_VERSION	"1.0"   /* Sphinx-3 file format version for mean/var */
#define NONE		-1
#define WORST_DIST	(int32)(0x80000000)

/** Subtract GMM component b (assumed to be positive) and saturate */
#ifdef FIXED_POINT
#define GMMSUB(a,b) \
	(((a)-(b) > a) ? (INT_MIN) : ((a)-(b)))
/** Add GMM component b (assumed to be positive) and saturate */
#define GMMADD(a,b) \
	(((a)+(b) < a) ? (INT_MAX) : ((a)+(b)))
#else
#define GMMSUB(a,b) ((a)-(b))
#define GMMADD(a,b) ((a)+(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif


#if defined(__STDC_VERSION__) && (__STDC_VERSION__ == 199901L)
#define LOGMATH_INLINE inline
#elif defined(__GNUC__)
#define LOGMATH_INLINE static inline
#elif defined(_MSC_VER)
#define LOGMATH_INLINE __inline
#else
#define LOGMATH_INLINE static
#endif

/* Allocate 0..159 for negated quantized mixture weights and 0..96 for
 * negated normalized acoustic scores, so that the combination of the
 * two (for a single mixture) can never exceed 255. */
#define MAX_NEG_MIXW 159 /**< Maximum negated mixture weight value. */
#define MAX_NEG_ASCR 96  /**< Maximum negated acoustic score value. */

/**
 * Quickly log-add two negated log probabilities.
 *
 * @param lmath The log-math object
 * @param mlx A negative log probability (0 < mlx < 255)
 * @param mly A negative log probability (0 < mly < 255)
 * @return -log(exp(-mlx)+exp(-mly))
 *
 * We can do some extra-fast log addition since we know that
 * mixw+ascr is always less than 256 and hence x-y is also always less
 * than 256.  This relies on some cooperation from logmath_t which
 * will never produce a logmath table smaller than 256 entries.
 *
 * Note that the parameters are *negated* log probabilities (and
 * hence, are positive numbers), as is the return value.  This is the
 * key to the "fastness" of this function.
 */
LOGMATH_INLINE int
fast_logmath_add(logmath_t *lmath, int mlx, int mly)
{
    logadd_t *t = LOGMATH_TABLE(lmath);
    int d, r;

    /* d must be positive, obviously. */
    if (mlx > mly) {
        d = (mlx - mly);
        r = mly;
    }
    else {
        d = (mly - mlx);
        r = mlx;
    }

    return r - (((uint8 *)t->table)[d]);
}

#endif /* __TIED_MGAU_COMMON_H__ */
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * tmat.c
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1997 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * $Log: tmat.c,v $
 * Revision 1.1.1.1  2006/05/23 18:45:01  dhuggins
 * re-importation
 *
 * Revision 1.4  2005/11/14 16:14:34  dhuggins
 * Use LOG() instead of logs3() for loading tmats, makes startup
 * ***much*** faster.
 *
 * Revision 1.3  2005/10/10 14:50:35  dhuggins
 * Deal properly with empty transition matrices.
 *
 * Revision 1.2  2005/09/30 15:01:23  dhuggins
 * More robust tmat reading - read the tmat in accordance with the fixed s2 topology
 *
 * Revision 1.1  2005/09/29 21:51:19  dhuggins
 * Add support for Sphinx3 tmat files.  Amazingly enough, it Just Works
 * (but it isn't terribly robust)
 *
 * Revision 1.6  2005/07/05 13:12:39  dhdfu
 * Add new arguments to logs3_init() in some tests, main_ep
 *
 * Revision 1.5  2005/06/21 19:23:35  arthchan2003
 * 1, Fixed doxygen documentation. 2, Added $ keyword.
 *
 * Revision 1.5  2005/05/03 04:09:09  archan
 * Implemented the heart of word copy search. For every ci-phone, every word end, a tree will be allocated to preserve its pathscore.  This is different from 3.5 or below, only the best score for a particular ci-phone, regardless of the word-ends will be preserved at every frame.  The graph propagation will not collect unused word tree at this point. srch_WST_propagate_wd_lv2 is also as the most stupid in the century.  But well, after all, everything needs a start.  I will then really get the results from the search and see how it looks.
 *
 * Revision 1.4  2005/04/21 23:50:26  archan
 * Some more refactoring on the how reporting of structures inside kbcore_t is done, it is now 50% nice. Also added class-based LM test case into test-decode.sh.in.  At this moment, everything in search mode 5 is already done.  It is time to test the idea whether the search can really be used.
 *
 * Revision 1.3  2005/03/30 01:22:47  archan
 * Fixed mistakes in last updates. Add
 *
 * 
 * 20.Apr.2001  RAH (rhoughton@mediasite.com, ricky.houghton@cs.cmu.edu)
 *              Added tmat_free to free allocated memory 
 *
 * 29-Feb-2000	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added tmat_chk_1skip(), and made tmat_chk_uppertri() public.
 * 
 * 10-Dec-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Made tmat_dump() public.
 * 
 * 11-Mar-97	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Started based on original S3 implementation.
 */

/* System headers. */
#include <string.h>

/* SphinxBase headers. */
#include <sphinxbase/logmath.h>
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/bio.h>

/* Local headers. */
#include "tmat.h"
#include "hmm.h"
#include "vector.h"

#define TMAT_PARAM_VERSION		"1.0"


/**
 * Checks that no transition matrix in the given object contains backward arcs.
 * @returns 0 if successful, -1 if check failed.
 */
static int32 tmat_chk_uppertri(tmat_t *tmat, logmath_t *lmath);


/**
 * Checks that transition matrix arcs in the given object skip over
 * at most 1 state.  
 * @returns 0 if successful, -1 if check failed.  
 */

static int32 tmat_chk_1skip(tmat_t *tmat, logmath_t *lmath);


void
tmat_dump(tmat_t * tmat, FILE * fp)
{
    int32 i, src, dst;

    for (i = 0; i < tmat->n_tmat; i++) {
        fprintf(fp, "TMAT %d = %d x %d\n", i, tmat->n_state,
                tmat->n_state + 1);
        for (src = 0; src < tmat->n_state; src++) {
            for (dst = 0; dst <= tmat->n_state; dst++)
                fprintf(fp, " %12d", tmat->tp[i][src][dst]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    fflush(fp);
}


/*
 * Check model tprob matrices that they conform to upper-triangular assumption;
 * i.e. no "backward" transitions allowed.
 */
int32
tmat_chk_uppertri(tmat_t * tmat, logmath_t *lmath)
{
    int32 i, src, dst;

    /* Check that each tmat is upper-triangular */
    for (i = 0; i < tmat->n_tmat; i++) {
        for (dst = 0; dst < tmat->n_state; dst++)
            for (src = dst + 1; src < tmat->n_state; src++)
                if (tmat->tp[i][src][dst] < 255) {
                    E_ERROR("tmat[%d][%d][%d] = %d\n",
                            i, src, dst, tmat->tp[i][src][dst]);
                    return -1;
                }
    }

    return 0;
}


int32
tmat_chk_1skip(tmat_t * tmat, logmath_t *lmath)
{
    int32 i, src, dst;

    for (i = 0; i < tmat->n_tmat; i++) {
        for (src = 0; src < tmat->n_state; src++)
            for (dst = src + 3; dst <= tmat->n_state; dst++)
                if (tmat->tp[i][src][dst] < 255) {
                    E_ERROR("tmat[%d][%d][%d] = %d\n",
                            i, src, dst, tmat->tp[i][src][dst]);
                    return -1;
                }
    }

    return 0;
}


tmat_t *
tmat_init(char const *file_name, logmath_t *lmath, float64 tpfloor, int32 breport)
{
    char tmp;
    int32 n_src, n_dst, n_tmat;
    FILE *fp;
    int32 byteswap, chksum_present;
    uint32 chksum;
    float32 **tp;
    int32 i, j, k, tp_per_tmat;
    char **argname, **argval;
    tmat_t *t;


    if (breport) {
        E_INFO("Reading HMM transition probability matrices: %s\n",
               file_name);
    }

    t = (tmat_t *) ckd_calloc(1, sizeof(tmat_t));

    if ((fp = fopen(file_name, "rb")) == NULL)
        E_FATAL_SYSTEM("Failed to open transition file '%s' for reading", file_name);

    /* Read header, including argument-value info and 32-bit byteorder magic */
    if (bio_readhdr(fp, &argname, &argval, &byteswap) < 0)
        E_FATAL("Failed to read header from file '%s'\n", file_name);

    /* Parse argument-value list */
    chksum_present = 0;
    for (i = 0; argname[i]; i++) {
        if (strcmp(argname[i], "version") == 0) {
            if (strcmp(argval[i], TMAT_PARAM_VERSION) != 0)
                E_WARN("Version mismatch(%s): %s, expecting %s\n",
                       file_name, argval[i], TMAT_PARAM_VERSION);
        }
        else if (strcmp(argname[i], "chksum0") == 0) {
            chksum_present = 1; /* Ignore the associated value */
        }
    }
    bio_hdrarg_free(argname, argval);
    argname = argval = NULL;

    chksum = 0;

    /* Read #tmat, #from-states, #to-states, arraysize */
    if ((bio_fread(&n_tmat, sizeof(int32), 1, fp, byteswap, &chksum)
         != 1)
        || (bio_fread(&n_src, sizeof(int32), 1, fp, byteswap, &chksum) !=
            1)
        || (bio_fread(&n_dst, sizeof(int32), 1, fp, byteswap, &chksum) !=
            1)
        || (bio_fread(&i, sizeof(int32), 1, fp, byteswap, &chksum) != 1)) {
        E_FATAL("Failed to read header from '%s'\n", file_name);
    }
    if (n_tmat >= MAX_INT16)
        E_FATAL("%s: Number of transition matrices (%d) exceeds limit (%d)\n", file_name,
                n_tmat, MAX_INT16);
    t->n_tmat = n_tmat;
    
    if (n_dst != n_src + 1)
        E_FATAL("%s: Unsupported transition matrix. Number of source states (%d) != number of target states (%d)-1\n", file_name,
                n_src, n_dst);
    t->n_state = n_src;

    if (i != t->n_tmat * n_src * n_dst) {
        E_FATAL
            ("%s: Invalid transitions. Number of coefficients (%d) doesn't match expected array dimension: %d x %d x %d\n",
             file_name, i, t->n_tmat, n_src, n_dst);
    }

    /* Allocate memory for tmat data */
    t->tp = ckd_calloc_3d(t->n_tmat, n_src, n_dst, sizeof(***t->tp));

    /* Temporary structure to read in the float data */
    tp = ckd_calloc_2d(n_src, n_dst, sizeof(**tp));

    /* Read transition matrices, normalize and floor them, and convert to log domain */
    tp_per_tmat = n_src * n_dst;
    for (i = 0; i < t->n_tmat; i++) {
        if (bio_fread(tp[0], sizeof(float32), tp_per_tmat, fp,
                      byteswap, &chksum) != tp_per_tmat) {
            E_FATAL("Failed to read transition matrix %d from '%s'\n", i, file_name);
        }

        /* Normalize and floor */
        for (j = 0; j < n_src; j++) {
            if (vector_sum_norm(tp[j], n_dst) == 0.0)
                E_WARN("Normalization failed for transition matrix %d from state %d\n",
                       i, j);
            vector_nz_floor(tp[j], n_dst, tpfloor);
            vector_sum_norm(tp[j], n_dst);

            /* Convert to logs3. */
            for (k = 0; k < n_dst; k++) {
                int ltp;
#if 0 /* No, don't do this!  It will subtly break 3-state HMMs. */
                /* For these ones, we floor them even if they are
                 * zero, otherwise HMM evaluation goes nuts. */
                if (k >= j && k-j < 3 && tp[j][k] == 0.0f)
                    tp[j][k] = tpfloor;
#endif
                /* Log and quantize them. */
                ltp = -logmath_log(lmath, tp[j][k]) >> SENSCR_SHIFT;
                if (ltp > 255) ltp = 255;
                t->tp[i][j][k] = (uint8)ltp;
            }
        }
    }

    ckd_free_2d(tp);

    if (chksum_present)
        bio_verify_chksum(fp, byteswap, chksum);

    if (fread(&tmp, 1, 1, fp) == 1)
        E_ERROR("Non-empty file beyond end of data\n");

    fclose(fp);

    if (tmat_chk_uppertri(t, lmath) < 0)
        E_FATAL("Tmat not upper triangular\n");
    if (tmat_chk_1skip(t, lmath) < 0)
        E_FATAL("Topology not Left-to-Right or Bakis\n");

    return t;
}

void
tmat_report(tmat_t * t)
{
    E_INFO_NOFN("Initialization of tmat_t, report:\n");
    E_INFO_NOFN("Read %d transition matrices of size %dx%d\n",
                t->n_tmat, t->n_state, t->n_state + 1);
    E_INFO_NOFN("\n");

}

/* 
 *  RAH, Free memory allocated in tmat_init ()
 */
void
tmat_free(tmat_t * t)
{
    if (t) {
        if (t->tp)
            ckd_free_3d(t->tp);
        ckd_free(t);
    }
}
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/*
 * tmat.h
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1997 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * $Log: tmat.h,v $
 * Revision 1.1.1.1  2006/05/23 18:45:03  dhuggins
 * re-importation
 *
 * Revision 1.1  2005/09/29 21:51:19  dhuggins
 * Add support for Sphinx3 tmat files.  Amazingly enough, it Just Works
 * (but it isn't terribly robust)
 *
 * Revision 1.9  2005/06/21 19:23:35  arthchan2003
 * 1, Fixed doxygen documentation. 2, Added $ keyword.
 *
 * Revision 1.6  2005/06/13 04:02:56  archan
 * Fixed most doxygen-style documentation under libs3decoder.
 *
 * Revision 1.5  2005/05/03 04:09:09  archan
 * Implemented the heart of word copy search. For every ci-phone, every word end, a tree will be allocated to preserve its pathscore.  This is different from 3.5 or below, only the best score for a particular ci-phone, regardless of the word-ends will be preserved at every frame.  The graph propagation will not collect unused word tree at this point. srch_WST_propagate_wd_lv2 is also as the most stupid in the century.  But well, after all, everything needs a start.  I will then really get the results from the search and see how it looks.
 *
 * Revision 1.4  2005/04/21 23:50:26  archan
 * Some more refactoring on the how reporting of structures inside kbcore_t is done, it is now 50% nice. Also added class-based LM test case into test-decode.sh.in.  At this moment, everything in search mode 5 is already done.  It is time to test the idea whether the search can really be used.
 *
 * Revision 1.3  2005/03/30 01:22:47  archan
 * Fixed mistakes in last updates. Add
 *
 * 
 * 20.Apr.2001  RAH (rhoughton@mediasite.com, ricky.houghton@cs.cmu.edu)
 *              Added tmat_free to free allocated memory 
 *
 * 29-Feb-2000	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added tmat_chk_1skip(), and made tmat_chk_uppertri() public.
 * 
 * 10-Dec-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added tmat_dump().
 * 
 * 11-Mar-97	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Started based on original S3 implementation.
 */


#ifndef _S3_TMAT_H_
#define _S3_TMAT_H_

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/logmath.h>

/** \file tmat.h
 *  \brief Transition matrix data structure.
 */
#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* Fool Emacs into not indenting things. */
#endif

/**
 * \struct tmat_t
 * \brief Transition matrix data structure.  All phone HMMs are assumed to have the same
 * topology.
 */
typedef struct {
    uint8 ***tp;	/**< The transition matrices; kept in the same scale as acoustic scores;
			   tp[tmatid][from-state][to-state] */
    int16 n_tmat;	/**< Number matrices */
    int16 n_state;	/**< Number source states in matrix (only the emitting states);
			   Number destination states = n_state+1, it includes the exit state */
} tmat_t;


/** Initialize transition matrix */

tmat_t *tmat_init (char const *tmatfile,/**< In: input file */
		   logmath_t *lmath,    /**< In: log math parameters */
		   float64 tpfloor,	/**< In: floor value for each non-zero transition probability */
		   int32 breport      /**< In: whether reporting the process of tmat_t  */
    );
					    


/** Dumping the transition matrix for debugging */

void tmat_dump (tmat_t *tmat,  /**< In: transition matrix */
		FILE *fp       /**< In: file pointer */
    );	


/**
 * RAH, add code to remove memory allocated by tmat_init
 */

void tmat_free (tmat_t *t /**< In: transition matrix */
    );

/**
 * Report the detail of the transition matrix structure. 
 */
void tmat_report(tmat_t *t /**< In: transition matrix*/
    );

#if 0
{ /* Stop indent from complaining */
#endif
#ifdef __cplusplus
}
#endif

#endif
/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * vector.c
 *
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1997 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * 
 * 22-Nov-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Imported from s3.2, for supporting s3 format continuous
 * 		acoustic models.
 * 
 * 10-Mar-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added vector_accum(), vector_vqlabel(), and vector_vqgen().
 * 
 * 09-Mar-1999	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added vector_is_zero(), vector_cmp(), and vector_dist_eucl().
 * 		Changed the name vector_dist_eval to vector_dist_maha.
 * 
 * 07-Oct-98	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Added distance computation related functions.
 * 
 * 12-Nov-95	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University.
 * 		Copied from Eric Thayer.
 */

/* System headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/* SphinxBase headers. */
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/bitvec.h>

/* Local headers. */
#include "vector.h"

#if (WIN32)
#define srandom	srand
#define random	rand
#endif


float64
vector_sum_norm(float32 * vec, int32 len)
{
    float64 sum, f;
    int32 i;

    sum = 0.0;
    for (i = 0; i < len; i++)
        sum += vec[i];

    if (sum != 0.0) {
        f = 1.0 / sum;
        for (i = 0; i < len; i++)
            vec[i] *= f;
    }

    return sum;
}


void
vector_floor(float32 * vec, int32 len, float64 flr)
{
    int32 i;

    for (i = 0; i < len; i++)
        if (vec[i] < flr)
            vec[i] = (float32) flr;
}


void
vector_nz_floor(float32 * vec, int32 len, float64 flr)
{
    int32 i;

    for (i = 0; i < len; i++)
        if ((vec[i] != 0.0) && (vec[i] < flr))
            vec[i] = (float32) flr;
}


void
vector_print(FILE * fp, vector_t v, int32 dim)
{
    int32 i;

    for (i = 0; i < dim; i++)
        fprintf(fp, " %11.4e", v[i]);
    fprintf(fp, "\n");
    fflush(fp);
}


int32
vector_is_zero(float32 * vec, int32 len)
{
    int32 i;

    for (i = 0; (i < len) && (vec[i] == 0.0); i++);
    return (i == len);          /* TRUE iff all mean values are 0.0 */
}
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/*
 * vector.h -- vector routines.
 * 
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1997 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 */


#ifndef __VECTOR_H__
#define __VECTOR_H__

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>

typedef float32 *vector_t;

/*
 * The reason for some of the "trivial" routines below is that they could be OPTIMIZED for SPEED
 * at some point.
 */


/* Floor all elements of v[0..dim-1] to min value of f */
void vector_floor(vector_t v, int32 dim, float64 f);


/* Floor all non-0 elements of v[0..dim-1] to min value of f */
void vector_nz_floor(vector_t v, int32 dim, float64 f);


/*
 * Normalize the elements of the given vector so that they sum to 1.0.  If the sum is 0.0
 * to begin with, the vector is left untouched.  Return value: The normalization factor.
 */
float64 vector_sum_norm(vector_t v, int32 dim);


/* Print vector in one line, in %11.4e format, terminated by newline */
void vector_print(FILE *fp, vector_t v, int32 dim);


/* Return TRUE iff given vector is all 0.0 */
int32 vector_is_zero (float32 *vec,	/* In: Vector to be checked */
		      int32 len);	/* In: Length of above vector */

#endif /* VECTOR_H */ 


/*
 * Log record.  Maintained by RCS.
 *
 * Revision 1.1.1.1  2004/03/29 20:29:40  rkm
 *
 *
 * Revision 1.1.1.1  2004/03/05 16:55:44  rkm
 *
 *
 * Revision 1.1.1.1  2000/02/28 18:05:54  rkm
 * Imported Sources
 *
 * Revision 1.1.1.1  1999/11/23 20:24:18  rkm
 * imported sources
 *
 * Revision 1.2  1995/10/09  20:55:35  eht
 * Changes for prim_type.h
 *
 * Revision 1.1  1995/08/15  13:44:14  eht
 * Initial revision
 *
 *
 */
