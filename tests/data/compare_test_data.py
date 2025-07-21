# Ellip is licensed under The 3-Clause BSD, see LICENSE.
# Copyright 2025 Sira Pornsiriprasert <code@psira.me>

"""
A script to compare two test data files

Take two files. Parse the files as CSV and extract the last value of each row.
Compare those values using symmetric relative error. Report descriptive statistics. 
"""

import pandas as pd
import numpy as np

boost_df = pd.read_csv("/home/psira/Code/rust/ellip/tests/data/boost/elliprf_data.csv", header=None)
wolfram_df = pd.read_csv("/home/psira/Code/rust/ellip/tests/data/wolfram/elliprf_data.csv", header=None)

def compare_dataframes_rerr(df1, df2):
    """
    Compare two dataframes using relative error.
    RERR = max(|a - b| / |a|, |a - b| / |b|)
    
    Returns statistics for each column: mean, median, p99, and max RERR.
    """
    
    # Ensure dataframes have same shape
    if df1.shape != df2.shape:
        raise ValueError(f"DataFrames must have same shape. Got {df1.shape} and {df2.shape}")
    
    # Ensure same columns
    if not df1.columns.equals(df2.columns):
        raise ValueError("DataFrames must have same column names")
    
    results = {}
    
    for col in df1.columns:
        # Get values for current column
        a = df1[col].values
        b = df2[col].values
        
        # Calculate relative error: max(|a-b|/|a|, |a-b|/|b|)
        abs_diff = np.abs(a - b)
        
        # Handle division by zero
        rerr1 = np.where(a == 0, np.inf, abs_diff / np.abs(a))
        rerr2 = np.where(b == 0, np.inf, abs_diff / np.abs(b))
        
        rerr = np.maximum(rerr1, rerr2)
        
        # Where both a and b are 0, set rerr to 0
        rerr = np.where((a == 0) & (b == 0), 0, rerr) / 2.22e-16

        rerr = np.maximum(0, rerr)
        
        # Remove any NaN or infinite values for statistics calculation
        rerr_clean = rerr[np.isfinite(rerr)]
        
        if len(rerr_clean) == 0:
            results[col] = {
                'mean': np.nan,
                'median': np.nan,
                'p99': np.nan,
                'max': np.nan,
                'count_valid': 0
            }
        else:
            results[col] = {
                'mean': np.mean(rerr_clean),
                'median': np.median(rerr_clean),
                'p99': np.percentile(rerr_clean, 99),
                'max': np.max(rerr_clean),
                'count_valid': len(rerr_clean)
            }
    
    return pd.DataFrame(results).T

print(compare_dataframes_rerr(boost_df, wolfram_df))