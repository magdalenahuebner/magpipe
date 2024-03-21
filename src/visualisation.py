import pandas as pd
import logging
import matplotlib.pyplot as plt


def plot_qc_samples(df, p_quant=0.6, save_path=None):
    """
    Plots the fractions of quantified phosphosites for each sample and either displays the plot or saves it to a file.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data with samples as columns.
    - p_quant (float): The quantification threshold, highlighted on the plot.
    - save_path (str, optional): Path to save the plot. If None, the plot is displayed directly.
    """
    # Calculate the fraction of quantified phosphosites for each sample
    df_quant = 1 - df.isnull().mean()
    
    # Create a bar plot
    plt.figure(figsize=(10, 6))
    df_quant.plot(kind='bar')
    plt.axhline(y=p_quant, color='r', linestyle='--')
    plt.title('Fraction of Quantified Phosphosites per Sample')
    plt.xlabel('Samples')
    plt.ylabel('Fraction Quantified')
    plt.xticks(rotation=90, fontsize=8)  # Rotate and resize sample names for better readability
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels

    # Save or show the plot based on save_path
    if save_path:
        plt.savefig(save_path)
        logging.info(f"Plot saved to {save_path}")
    else:
        plt.show()


def plot_normalisation_boxplots(df, df_norm, save_path=None):
    """
    Generates boxplots for each sample in the DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data. Rows represent features, and columns represent samples.
    - title (str): Title of the plot.
    """
    # Create boxplots for the original data
    plt.figure(figsize=(12, 6))
    df.boxplot(rot=90, fontsize=8)
    plt.title('Before Normalisation')
    plt.ylabel('Area under the Peak (AUP)')  # Rotate and resize sample names for better readability
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels

    # Create boxplots for the normalised data
    plt.figure(figsize=(12, 6))
    df_norm.boxplot(rot=90, fontsize=8)
    plt.title('After Normalisation')
    plt.ylabel('Area under the Peak (AUP)')  # Rotate and resize sample names for better readability
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels

    if save_path:
        plt.savefig(save_path)
        logging.info(f"QC plots saved to {save_path}")
    else:
        plt.show()
