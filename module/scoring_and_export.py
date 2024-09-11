# NMFF-AFM: normal mode flexible fitting of proteins conformations to AFM images
#
# Copyright (c) 2024 TamaLab
# Authors: Xuan Wu, Nagoya University
#          Osamu Miyashita, RIKEN
#          Florence Tama, Nagoya University/RIKEN
#
# References:
# Modeling Conformational Transitions of Biomolecules from Atomic Force
# Microscopy Images using Normal Mode Analysis
# Xuan Wu, Osamu Miyashita, and Florence Tama
# https://doi.org/10.1021/acs.jpcb.4c04189
#
#
import openpyxl
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def read_log(log_path, rmsd_calculation):
    """
    Reads a log file and processes the data based on whether RMSD calculation is required
    :param log_path: The path to the log file (Excel format) to be read
    :param rmsd_calculation: A flag indicating whether RMSD calculation is enabled.
        If "yes", the function reads additional column related to RMSD-Reference
    :return: A DataFrame containing the processed log data
    """
    if rmsd_calculation == "yes":
        read_sheet = pd.read_excel(log_path, usecols="A:G")
        read_sheet["step"] = range(len(read_sheet.index))
        read_sheet = read_sheet.rename(columns={"CC of this iter": "cc"})
        read_sheet["cc change"] = read_sheet["cc"].diff()
    else:
        read_sheet = pd.read_excel(log_path, usecols="A:F")
        read_sheet["step"] = range(len(read_sheet.index))
        read_sheet = read_sheet.rename(columns={"CC of this iter": "cc"})
        read_sheet["cc change"] = read_sheet["cc"].diff()

    return read_sheet


def save_time_series_plot(read_sheet, time_series_fig_path, time_series_annotated_fig_path):
    """
    Generates and saves time series plots for RMSD, correlation coefficient,
    and mode data from a simulation run with the data related to RMSD-Reference.
    Highlights significant changes and performs an exponential fit to the data
    :param read_sheet: A DataFrame containing the simulation data, with the column of RMSD-Reference
    :param time_series_fig_path: The file path where the main figure with time series will be saved
    :param time_series_annotated_fig_path: The file path where the annotated plot will be saved
    :return: The number of steps at which the correlation coefficient change decays to specific factors,
            primarily used to determine when certain conditions are met within the simulation.
    """
    global num_of_step
    fig, ax = plt.subplots(4, figsize=(6, 12))
    A_RMSD = read_sheet.loc[:, ["RMSD to Initial", "RMSD to Reference"]]

    sns.lineplot(data=A_RMSD, ax=ax[0])
    ax[0].set_ylabel('rmsd')
    ax[0].set_xlabel('')

    sns.lineplot(read_sheet, x="step", y="cc", ax=ax[1])
    ax[1].set_xlabel('')

    sns.lineplot(read_sheet, x="step", y="cc change", ax=ax[2])
    ax[2].axhline(y=0, linestyle=':')
    ax[2].set_xlabel('')

    sns.lineplot(read_sheet, x="step", y="Largest mode", ax=ax[3])
    fig.align_ylabels()
    fig.tight_layout()
    plt.savefig(time_series_fig_path)

    # Check if "RMSD to Reference" column exists and DataFrame is not empty
    if "RMSD to Reference" in A_RMSD.columns and not A_RMSD.empty:
        # Plots with annotations
        min_index = A_RMSD["RMSD to Reference"].idxmin()

        # Ensure the index is valid
        if isinstance(min_index, int):
            x_min = A_RMSD.index[min_index]  # Assuming your index is numeric or datetime
            y_min = A_RMSD["RMSD to Reference"].iloc[min_index]

            # Plot the line
            plt.plot([x_min, x_min], [y_min, y_min], 'r--')
        else:
            print("Error: Unable to find minimum index of RMSD to Reference, maybe reference.pdb does not exist.")
    else:
        print("Warning: Column 'RMSD to Reference' not found or DataFrame is empty. Skipping line plot.")

    # Find the first index where 'cc change' is negative
    first_negative_index = (read_sheet["cc change"] < 0).idxmax()

    # Get the coordinates of the first negative point
    x_value = read_sheet.loc[first_negative_index, "step"]
    y_value = read_sheet.loc[first_negative_index, "cc change"]
    txt = f"({x_value}, {y_value:.{2}g})"
    ax[2].scatter(x_value, y_value, color="red", marker="v", label=f"1st neg {txt}")
    y_value = read_sheet.loc[x_value, 'RMSD to Reference']
    ax[0].scatter(x_value, y_value, color="red", marker="v", label=f"1st neg {txt}")

    # Define the exponential function
    def exponential_func(x, a, b):
        return a * np.exp(b * (x - 1))

    # Fit the cc change data to the exponential function
    x = read_sheet["step"].iloc[1:]  # the first value is None
    y = read_sheet["cc change"].iloc[1:]
    params, covariance = curve_fit(exponential_func, x, y, p0=[1, -1])
    print('exp fit: ', params, covariance)
    a_fit, b_fit = params
    ax[2].plot(
        x,
        exponential_func(x, a_fit, b_fit),
        label=f"y={a_fit:.{2}g}exp({b_fit:.{2}g}(x-1))",
    )

    # find the point at the specific decays
    p_list = [0.05, 0.03, 0.01]
    for p in p_list:
        u = int(np.ceil(np.log(p) / b_fit + 1))
        m = A_RMSD.index.max()
        if u > m:
            print(f'Warning: Diff C is expected to decay by the factor {p} at step {u}, but this run stopped at {m}')
            u = m
        if p == 0.03:
            num_of_step = u
        v = exponential_func(u, a_fit, b_fit)
        txt = f"{p}: ({u}, {v:.{1}e})"
        ax[2].scatter(u, v, marker="o", label=txt)

        w = A_RMSD.loc[u, 'RMSD to Reference']
        txt = f"{p}: ({u}, {w:.{3}g})"
        ax[0].scatter(u, w, marker='o', label=txt)
    ax[2].legend(loc='center left', bbox_to_anchor=(1., 0.5))
    ax[0].legend(loc='center left', bbox_to_anchor=(1., 0.5))

    plt.savefig(time_series_annotated_fig_path, bbox_inches='tight')
    return num_of_step


def save_time_series_plot_no_reference(read_sheet, time_series_fig_path, time_series_annotated_fig_path):
    """
    Generates and saves time series plots for RMSD, correlation coefficient,
    and mode data from a simulation run without the data related to RMSD-Reference.
    Highlights significant changes and performs an exponential fit to the data
    :param read_sheet: A DataFrame containing the simulation data, without the column of RMSD-Reference
    :param time_series_fig_path: The file path where the main figure with time series will be saved
    :param time_series_annotated_fig_path: The file path where the annotated plot will be saved
    :return: The number of steps at which the correlation coefficient change decays to specific factors,
            primarily used to determine when certain conditions are met within the simulation.
    """
    global num_of_step
    fig, ax = plt.subplots(4, figsize=(6, 12))
    A_RMSD = read_sheet.loc[:, ["RMSD to Initial"]]

    sns.lineplot(data=A_RMSD, ax=ax[0])
    ax[0].set_ylabel('rmsd')
    ax[0].set_xlabel('')

    sns.lineplot(read_sheet, x="step", y="cc", ax=ax[1])
    ax[1].set_xlabel('')

    sns.lineplot(read_sheet, x="step", y="cc change", ax=ax[2])
    ax[2].axhline(y=0, linestyle=':')
    ax[2].set_xlabel('')

    sns.lineplot(read_sheet, x="step", y="Largest mode", ax=ax[3])
    fig.align_ylabels()
    fig.tight_layout()
    plt.savefig(time_series_fig_path)

    # Plots with annotations
    min_index = A_RMSD["RMSD to Initial"].idxmin()

    # Ensure the index is valid
    if isinstance(min_index, int):
        x_min = A_RMSD.index[min_index]  # Assuming your index is numeric or datetime
        y_min = A_RMSD["RMSD to Initial"].iloc[min_index]

        # Plot the line
        plt.plot([x_min, x_min], [y_min, y_min], 'r--')

    # Find the first index where 'cc change' is negative
    first_negative_index = (read_sheet["cc change"] < 0).idxmax()

    # Get the coordinates of the first negative point
    x_value = read_sheet.loc[first_negative_index, "step"]
    y_value = read_sheet.loc[first_negative_index, "cc change"]
    txt = f"({x_value}, {y_value:.{2}g})"
    ax[2].scatter(x_value, y_value, color="red", marker="v", label=f"1st neg {txt}")

    # Define the exponential function
    def exponential_func(x, a, b):
        return a * np.exp(b * (x - 1))

    # Fit the cc change data to the exponential function
    x = read_sheet["step"].iloc[1:]  # the first value is None
    y = read_sheet["cc change"].iloc[1:]
    params, covariance = curve_fit(exponential_func, x, y, p0=[1, -1])
    print('exp fit: ', params, covariance)
    a_fit, b_fit = params
    ax[2].plot(
        x,
        exponential_func(x, a_fit, b_fit),
        label=f"y={a_fit:.{2}g}exp({b_fit:.{2}g}(x-1))",
    )

    # find the point at the specific decays
    p_list = [0.05, 0.03, 0.01]
    for p in p_list:
        u = int(np.ceil(np.log(p) / b_fit + 1))
        m = A_RMSD.index.max()
        if u > m:
            print(f'Warning: Diff C is expected to decay by the factor {p} at step {u}, but this run stopped at {m}')
            u = m
        if p == 0.03:
            num_of_step = u
        v = exponential_func(u, a_fit, b_fit)
        txt = f"{p}: ({u}, {v:.{1}e})"
        ax[2].scatter(u, v, marker="o", label=txt)

        w = A_RMSD.loc[u, 'RMSD to Initial']
        txt = f"{p}: ({u}, {w:.{3}g})"
        ax[0].scatter(u, w, marker='o', label=txt)
    ax[2].legend(loc='center left', bbox_to_anchor=(1., 0.5))
    ax[0].legend(loc='center left', bbox_to_anchor=(1., 0.5))

    plt.savefig(time_series_annotated_fig_path, bbox_inches='tight')
    return num_of_step


def plot_and_save_gain_figure(read_sheet, step_pick, result_fig_path):
    """
    Plots time series data from a given DataFrame and saves the resulting figure, with the column related to RMSD-Reference
    :param read_sheet: A DataFrame containing the data to be plotted, including columns for "RMSD to Initial",
                    "RMSD to Reference", "cc" (correlation coefficient), and "Largest mode".
    :param step_pick: The step index at which specific data points will be highlighted on the plots
                    with scatter markers
    :param result_fig_path: The path where the generated figure will be saved
    :return: None
    """
    fig, ax = plt.subplots(3, figsize=(6, 10))
    A_RMSD = read_sheet.loc[:, ["RMSD to Initial", "RMSD to Reference"]]
    sns.lineplot(data=A_RMSD, ax=ax[0])

    ax[0].scatter(
        x=step_pick, y=read_sheet.loc[step_pick, "RMSD to Initial"], marker="o"
    )
    ax[0].scatter(
        x=step_pick, y=read_sheet.loc[step_pick, "RMSD to Reference"], marker="o"
    )
    sns.lineplot(read_sheet, x="step", y="cc", ax=ax[1])
    ax[1].scatter(
        x=step_pick, y=read_sheet.loc[step_pick, "cc"], marker="o"
    )
    sns.lineplot(read_sheet, x="step", y="Largest mode", ax=ax[2])
    ax[2].scatter(x=step_pick, y=read_sheet.loc[step_pick, "Largest mode"], marker="o")

    plt.savefig(result_fig_path)


def plot_and_save_gain_figure_no_reference(read_sheet, step_pick, result_fig_path):
    """
    Plots time series data from a given DataFrame and saves the resulting figure, without the column related to RMSD-Reference
    :param read_sheet: A DataFrame containing the data to be plotted, including columns for "RMSD to Initial",
                    "cc" (correlation coefficient), and "Largest mode".
    :param step_pick: The step index at which specific data points will be highlighted on the plots
                    with scatter markers
    :param result_fig_path: The path where the generated figure will be saved
    :return: None
    """
    fig, ax = plt.subplots(3, figsize=(6, 10))
    A_RMSD = read_sheet.loc[:, ["RMSD to Initial"]]
    sns.lineplot(data=A_RMSD, ax=ax[0])

    ax[0].scatter(
        x=step_pick, y=read_sheet.loc[step_pick, "RMSD to Initial"], marker="o"
    )
    sns.lineplot(read_sheet, x="step", y="cc", ax=ax[1])
    ax[1].scatter(
        x=step_pick, y=read_sheet.loc[step_pick, "cc"], marker="o"
    )
    sns.lineplot(read_sheet, x="step", y="Largest mode", ax=ax[2])
    ax[2].scatter(x=step_pick, y=read_sheet.loc[step_pick, "Largest mode"], marker="o")

    plt.savefig(result_fig_path)


def write_result_step_to_xlsx(io, result):
    """
    Writes the result of a specific step to an Excel file, updating a cell
    based on the result value
    :param io: The path to the Excel file or a file-like object where the result will be written
    :param result: The step result value to be written. This value is used to determine the row position
                in the sheet and is formatted as a string prefixed with 's' before being written.
    :return: None
    """
    full_result = f"s{result}"
    wb = openpyxl.load_workbook(io)
    sheet = wb.active
    cell = sheet.cell(row=int(result + 2), column=9)
    cell.value = full_result

    wb.save(io)
    print("Result written to log file.")


def scoring_conformations(work_folder, log_path, calculate_rmsd_reference):
    """
    Scores conformations based on time series data, generates plots, and saves results to specified path
    :param work_folder: The directory path where output figures will be saved
    :param log_path: The path to the log file that contains the time series data
    :param calculate_rmsd_reference: A flag ("yes" or "no") indicating whether to calculate and use RMSD to a reference conformation
                                    - If "yes", the function generates figures and calculations that include RMSD-Reference data
                                    - If "no", figures and calculations are generated without RMSD-Reference data
    :return: The index of the turning point in the time series data
    """
    read_sheet = read_log(log_path, calculate_rmsd_reference)

    time_series_fig = Path(work_folder) / "time_series.png"
    result_fig = Path(work_folder) / "result_figure.png"
    time_series_annotated_fig = Path(work_folder) / "time_series_annotated.png"

    if calculate_rmsd_reference == "yes":
        turning_point = save_time_series_plot(read_sheet, time_series_fig, time_series_annotated_fig)
        plot_and_save_gain_figure(read_sheet, turning_point, result_fig)
        write_result_step_to_xlsx(log_path, turning_point)
    else:
        turning_point = save_time_series_plot_no_reference(read_sheet, time_series_fig, time_series_annotated_fig)
        plot_and_save_gain_figure_no_reference(read_sheet, turning_point, result_fig)
        write_result_step_to_xlsx(log_path, turning_point)

    print("All the figures generated and save.")
    return turning_point
