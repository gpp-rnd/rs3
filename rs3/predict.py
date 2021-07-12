# AUTOGENERATED! DO NOT EDIT! File to edit: 04_predict.ipynb (unless otherwise specified).

__all__ = ['predict_seq_tracr', 'combine_target_seq_scores', 'predict']

# Cell
import pandas as pd
import warnings

from .seq import predict_seq
from .targetdata import (build_translation_overlap_df,
                            build_transcript_aa_seq_df)
from .predicttarg import predict_target

# Cell
def predict_seq_tracr(design_df, tracr, context_col, ref_tracrs):
    if not tracr in ref_tracrs:
        raise ValueError('tracrRNA must be one of ' + ','.join(ref_tracrs))
    design_df['RS3 Sequence Score (' + tracr + ' tracr)'] = predict_seq(design_df[context_col], sequence_tracr=tracr)

def combine_target_seq_scores(design_df, tracr):
    design_df['RS3 Sequence (' + tracr + ' tracr) + Target Score'] = \
        design_df['RS3 Sequence Score (' + tracr + ' tracr)'] + \
        design_df['RS3 Target Score']

def predict(design_df, tracr=None, target=False,
            aa_seq_file=None, domain_file=None,
            target_id_cols=None,
            context_col='sgRNA Context Sequence',
            transcript_id_col='Target Transcript',
            transcript_base_col='Transcript Base',
            transcript_len_col='Target Total Length',
            n_jobs=1):
    """Make predictions using RS3

    :param design_df: DataFrame
    :param tracr: str or list
    :param target: bool, whether to include target scores
    :param aa_seq_file: str, path to precomputed amino acid sequences
    :param domain_file: str, path to precomputed domain file
    :param target_id_cols: list or None
    :param context_col: str
    :param transcript_id_col: str
    :param transcript_base_col: str
    :param transcript_len_col: str
    :param n_jobs: int
    :return: DataFram
    """
    out_df = design_df.copy()
    ref_tracrs = ['Hsu2013', 'Chen2013']
    if type(tracr) is str:
        predict_seq_tracr(out_df, tracr, context_col, ref_tracrs)
    else: # list
        for t in tracr:
            predict_seq_tracr(out_df, t, context_col, ref_tracrs)
    if target:
        out_df[transcript_base_col] = out_df[transcript_id_col].str.split('.', expand=True)[0]
        transcript_bases = pd.Series(out_df[transcript_base_col].unique())
        if aa_seq_file is None:
            aa_seq_df = build_transcript_aa_seq_df(out_df,
                                                   transcript_id_col=transcript_id_col,
                                                   transcript_len_col=transcript_len_col,
                                                   n_jobs=n_jobs)
        else:
            aa_seq_df = pd.read_parquet(aa_seq_file, engine='pyarrow',
                                        filters=[[(transcript_base_col, 'in', transcript_bases)]])
        missing_transcripts_aa = transcript_bases[~transcript_bases.isin(aa_seq_df[transcript_base_col])]
        if len(missing_transcripts_aa) > 0:
            warnings.warn('Missing amino acid sequences for transcripts: ' +
                          ','.join(missing_transcripts_aa))
            out_df['Missing translation information'] = out_df[transcript_base_col].isin(missing_transcripts_aa)
        if domain_file is None:
            domain_df = build_translation_overlap_df(aa_seq_df['id'].unique(), n_jobs=n_jobs)
        else:
            domain_df = pd.read_parquet(domain_file, engine='pyarrow',
                                        filters=[[(transcript_base_col, 'in', out_df[transcript_base_col].unique())]])
        # No warning for domain, since some transcripts aren't annotated with any domains
        out_df['RS3 Target Score'] = predict_target(design_df=out_df, aa_seq_df=aa_seq_df,
                                                    protein_domain_df=domain_df,
                                                    id_cols=target_id_cols)
        if type(tracr) is str:
            combine_target_seq_scores(out_df, tracr)
        else: # list
            for t in tracr:
                combine_target_seq_scores(out_df, t)
    return out_df