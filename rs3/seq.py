# AUTOGENERATED! DO NOT EDIT! File to edit: 00_seq.ipynb (unless otherwise specified).

__all__ = ['load_seq_model', 'featurize_context', 'predict_seq']

# Cell
import joblib
import sglearn
import pandas as pd
import os

# Cell
def load_seq_model():
    """Load rule set 3 sequence model"""
    model = joblib.load(os.path.join(os.path.dirname(__file__), 'RuleSet3.pkl'))
    return model

# Cell
def featurize_context(context_sequences, sequence_tracr='Hsu2013', ref_tracrs=None,
                      n_jobs=1):
    """Featurize context sequences

    :param context_sequences: list-like
    :param sequence_tracr: list-like or str
    :return: DataFrame, feature matrix
    """
    if ref_tracrs is None:
        ref_tracrs = ['Hsu2013', 'Chen2013']
    context_series = pd.Series(context_sequences)
    if not (context_series.str.len() == 30).all():
        raise  ValueError('All context sequences must be 30 nucleotides')
    featurized_sgrnas = sglearn.featurize_guides(context_sequences,
                                                 n_jobs=n_jobs)
    for tracr in ref_tracrs:
        if type(sequence_tracr) is str:
            featurized_sgrnas[tracr + ' tracr'] = int(sequence_tracr == tracr)
        else: # list-like
            featurized_sgrnas[tracr + ' tracr'] = ((pd.Series(sequence_tracr) == tracr)
                                                   .astype(int)
                                                   .to_list())
    return featurized_sgrnas

# Cell
def predict_seq(context_sequences, sequence_tracr='Hsu2013', ref_tracrs=None, n_jobs=1):
    """Predict the activity of context sequence for SpCas9 Knockout using sequence information only

    :param context_sequences: list of str
    :return: list of float, predictions
    """
    model = load_seq_model()
    print('Calculating sequence-based features')
    featurized_sgrnas = featurize_context(context_sequences, sequence_tracr=sequence_tracr, ref_tracrs=ref_tracrs,
                                          n_jobs=n_jobs)
    seq_predictions = model.predict(featurized_sgrnas)
    return seq_predictions