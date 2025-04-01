import numpy as np
import pandas as pd
from pathlib import Path
from taxontabletools2.utilities import simple_taxontable
from taxontabletools2.utilities import filter_taxontable
from taxontabletools2.utilities import taxon_richness_df
import plotly.graph_objects as go
import statsmodels.api as sm

def time_series_richness_with_ci(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df,
                                 selected_metadata, traits_df,
                                 selected_traits, user_settings, tool_settings):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    ## collect tool-specific settings
    selected_metadata = tool_settings['selected_metadata']
    taxonomic_level_1 = tool_settings['taxonomic_level_1']
    taxonomic_level_2 = tool_settings['taxonomic_level_2']
    n_bootstraps = tool_settings['n_bootstraps']
    ci_level = tool_settings['ci_level']

    richness_df = taxon_richness_df(taxon_table_df, samples, taxonomic_level_1)

    ## collect the number of categories
    samples_sorted = metadata_df['Sample'].values.tolist()

    # Prepare x-values and y-values (using actual sample names as x-values)
    x_values = samples_sorted  # Use the sample names instead of index numbers
    y_values = richness_df[samples_sorted].values.tolist()[0]

    # Plot original data points
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_values, y=y_values, name='Original Data', mode='markers', marker_color=user_settings['color_1'], marker=dict(size=user_settings['scatter_size']) ))

    # Calculate LOESS (Locally Estimated Scatterplot Smoothing) for the original data
    # Since x_values are now sample names, we need numerical indices for LOESS smoothing
    x_numeric = list(range(len(x_values)))
    loess_smoothed = sm.nonparametric.lowess(endog=y_values, exog=x_numeric, frac=0.3)
    loess_x_numeric = loess_smoothed[:, 0]
    loess_y = loess_smoothed[:, 1]

    # Map numerical x-values back to sample names
    loess_x = [x_values[int(i)] for i in loess_x_numeric]

    # Add LOESS-smoothed line to the plot
    fig.add_trace(go.Scatter(x=loess_x, y=loess_y, name='LOESS Smoothing', mode='lines', marker_color=user_settings['color_2'], marker=dict(size=user_settings['scatter_size']) ))

    # Bootstrapping to calculate confidence intervals
    bootstrapped_fits = []

    for _ in range(n_bootstraps):
        # Resample the data with replacement
        resample_indices = np.random.choice(range(len(x_numeric)), size=len(x_numeric), replace=True)
        resample_x_numeric = np.array(x_numeric)[resample_indices]
        resample_y = np.array(y_values)[resample_indices]

        # Fit LOESS on resampled data
        loess_bootstrap = sm.nonparametric.lowess(endog=resample_y, exog=resample_x_numeric, frac=0.3)

        # Interpolate to original x-numeric values
        interp_bootstrap_y = np.interp(x_numeric, loess_bootstrap[:, 0], loess_bootstrap[:, 1])
        bootstrapped_fits.append(interp_bootstrap_y)

    # Convert list to numpy array for easier calculations
    bootstrapped_fits = np.array(bootstrapped_fits)

    # Calculate the confidence intervals (CI) at each x-value
    lower_bound = np.percentile(bootstrapped_fits, (100 - ci_level) / 2, axis=0)
    upper_bound = np.percentile(bootstrapped_fits, 100 - (100 - ci_level) / 2, axis=0)

    # Add shaded region for the confidence interval
    fig.add_trace(go.Scatter(
        x=x_values + x_values[::-1],  # x values forward and backward for the fill area
        y=np.concatenate([upper_bound, lower_bound[::-1]]),  # upper bound followed by reversed lower bound
        fill='toself',
        fillcolor='rgba(0,100,80,0.2)',  # Semi-transparent color
        line=dict(color='rgba(255,255,255,0)'),
        showlegend=False,
        name=f'{ci_level}% Confidence Interval'
    ))

    # Update layout
    fig.update_xaxes(dtick='linear')
    fig.update_yaxes(rangemode='tozero')
    fig.update_layout(
        title=f'Richness ({taxonomic_level_1}) with LOESS Smoothing and {ci_level}% Confidence Interval',
        template=user_settings['template'],
        font_size=user_settings['font_size'],
        yaxis_title=f'Richness ({taxonomic_level_1})',
        yaxis_title_font=dict(size=user_settings['font_size']),
        yaxis_tickfont=dict(size=user_settings['font_size']),
        xaxis_title='Sample',
        xaxis_tickvals=x_values,  # Ensure sample names are shown as tick labels
        xaxis_ticktext=x_values,  # Use sample names as tick labels
        xaxis_tickfont=dict(size=user_settings['font_size'])
    )

    # if user_settings['show_yaxis'] == False:
    #    fig.update_xaxes(showticklabels=False)

    return fig



