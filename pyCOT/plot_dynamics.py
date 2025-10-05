from pyvis.network import Network        # Enables creating and visualizing interactive networks in HTML.
import networkx as nx                    # Library for creating, manipulating, and analyzing graphs.
import pandas as pd                      # Efficient handling and analysis of tabular data.
import matplotlib.pyplot as plt          # For generating static plots and visualizations.
import matplotlib.animation as animation # For creating animations using Matplotlib.
from matplotlib.colors import to_hex     # Converts colors to hexadecimal format (#RRGGBB).
import mplcursors                        # Adds interactive cursors to Matplotlib plots.
import webbrowser                        # Opens URLs or local files in the system's default web browser.
import os                                # Handles file and directory operations in the operating system.
from collections import defaultdict      # Dictionary that provides default values for missing keys.
import plotly.graph_objects as go        # Generates interactive plots and advanced visualizations using Plotly.
from itertools import combinations       # Generates all possible combinations of elements in an iterable.
import numpy as np                       # Library for numerical computing in Python.

from networkx.drawing.nx_agraph import graphviz_layout
import itertools 
from mpl_toolkits.mplot3d import Axes3D
 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import sympy as sp 
from scipy.optimize import linprog  
from collections import Counter  
import warnings 

import warnings
warnings.filterwarnings('ignore')

import sys                               # Provides access to system-specific parameters and functions.
sys.stdout.reconfigure(encoding='utf-8') # Reconfigures the standard output to use UTF-8 encoding, ensuring proper handling of special characters.
import tempfile                          # Provides utilities for creating temporary files and directories.

######################################################################################
# Plots the time series of ODE concentrations and abstractions
######################################################################################
import matplotlib.pyplot as plt
import os

def plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", 
                   title="Time Series of Concentrations", filename="time_series_plot.png",
                   show_grid=True, save_figure=True,
                   ax=None, show_fig=False):
    if 'Time' not in time_series.columns:
        raise ValueError("The DataFrame must include a 'Time' column for time values.")

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.get_figure()

    for species in time_series.columns:
        if species != 'Time':
            ax.plot(time_series['Time'], time_series[species], label=species)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(show_grid)
    ax.legend()

    if save_figure:
        os.makedirs("visualizations/plot_series_ode", exist_ok=True)
        filepath = os.path.join("visualizations", "plot_series_ode", filename)
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filepath}")
    if show_fig:
        plt.show()
    return fig, ax

######################################################################################## 
# Reaction-Diffusion Dynamics Plotting
########################################################################################
# Function to plot time series for each species in each grid cell
import matplotlib.pyplot as plt
import os

def plot_diffusion_time_series_2D(time, concentration_data, grid_shape, colors=None, 
                                 xlabel="Time", ylabel="Concentration", 
                                 main_title="Time Evolution of Concentration Profiles",
                                 legend_title="Species", cell_prefix='G',
                                 save_figure=True, filename="diffusion_time_series_2D.png"):
    """
    Plot 2D time series of concentration data from diffusion simulation and optionally save the figure.
    
    Parameters:
    -----------
    time : array-like
        Time vector for the x-axis
    concentration_data : dict
        Dictionary with species names as keys and 3D arrays (time, row, col) as values
    grid_shape : tuple (rows, cols)
        Shape of the spatial grid
    colors : dict, optional
        Color mapping for each species
    xlabel, ylabel : str
        Axis labels
    main_title : str
        Main title for the entire figure
    legend_title : str
        Title for the legend
    cell_prefix : str
        Prefix for cell labels (e.g., 'G' for G00, G01, etc.)
    save_figure : bool, optional
        Whether to save the figure. Default is True.
    filename : str, optional
        Name of the file to save. Default is "diffusion_time_series_2D.png".
    """
    
    rows, cols = grid_shape
    total_cells = rows * cols
    
    # Create figure with adjusted width for the legend
    fig, axs = plt.subplots(rows, cols, figsize=(15, 8), sharex=True, squeeze=False)
    axs_flat = axs.flatten()
    
    # Get default color cycle
    prop_cycle = plt.rcParams['axes.prop_cycle']
    default_colors = prop_cycle.by_key()['color']
    
    # Create color mapping
    species_list = list(concentration_data.keys())
    if colors is None:
        colors = {}
    
    color_mapping = {}
    for i, species in enumerate(species_list):
        color_mapping[species] = colors.get(species, default_colors[i % len(default_colors)])
    
    # Dictionary to store legend handles (one per species)
    legend_handles = {}
    
    # Plot data for each cell
    for cell_idx in range(total_cells):
        row = cell_idx // cols
        col = cell_idx % cols
        ax = axs_flat[cell_idx]
        
        # Plot each species in this cell
        for species in species_list:
            line = ax.plot(time, concentration_data[species][:, row, col], 
                          color=color_mapping[species], linewidth=2, label=species)
            
            # Store handle for legend (first occurrence only)
            if species not in legend_handles:
                legend_handles[species] = line[0]
        
        # Set cell title and labels
        cell_label = f'{cell_prefix}{row}{col}'
        ax.set_title(cell_label, fontsize=12, fontweight='bold')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
    
    # Create legend
    handles = [legend_handles[sp] for sp in species_list]
    legend = fig.legend(handles, species_list, 
                       loc='center right',
                       bbox_to_anchor=(0.98, 0.5),
                       fontsize=10, 
                       title=legend_title, 
                       title_fontsize=12)
    
    # Set main title
    plt.suptitle(main_title, fontsize=14, fontweight='bold')
    
    # Adjust layout to accommodate legend
    plt.tight_layout(rect=[0, 0, 0.88, 0.95])
    
    # Save the figure if requested
    if save_figure:
        # Create the directory if it doesn't exist
        os.makedirs("visualizations/plot_diffusion_time_series_2D", exist_ok=True)
        
        # Construct the full file path
        filepath = os.path.join("visualizations", "plot_diffusion_time_series_2D", filename)
        
        # Save the figure with high quality
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filepath}")
    
    plt.show()

# Function to plot heatmaps for a given species at selected time points
import numpy as np
import matplotlib.pyplot as plt
import mplcursors

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import mplcursors
import os

def plot_heatmaps_all_species_2D(t, X, species_names=None, time_indices=None, 
                                main_title="Evolution of Concentration Profiles", 
                                cmap='viridis', figsize_multiplier=3, bar_label="Concentration",
                                save_figure=True, filename="heatmaps_all_species_2D.png"):
    """
    Plot heatmaps for all species at multiple time points in a single figure and optionally save the figure.
    
    Parameters:
    -----------
    t : array-like
        Time vector
    X : dict
        Dictionary with species names as keys and 3D arrays (time, rows, cols) as values
    species_names : list, optional
        List of species to plot. If None, plots all species in X.
    time_indices : int or list of int, optional
        Single time index or list of time indices to display. 
        If None, shows 5 evenly spaced time points.
        If integer, shows that many evenly spaced time points.
    main_title : str
        Main title for the entire figure
    cmap : str
        Colormap to use for heatmaps
    figsize_multiplier : float
        Multiplier for figure size calculation
    bar_label : str
        Label for the color bar
    save_figure : bool, optional
        Whether to save the figure. Default is True.
    filename : str, optional
        Name of the file to save. Default is "heatmaps_all_species_2D.png".
    """
    
    if species_names is None:
        species_names = list(X.keys())
    
    n_species = len(species_names)
    
    # Manejo de índices de tiempo
    if time_indices is None:
        time_indices = np.linspace(0, len(t)-1, 5, dtype=int)
    elif isinstance(time_indices, int):
        time_indices = np.linspace(0, len(t)-1, time_indices, dtype=int)
    elif not hasattr(time_indices, '__len__'):
        raise ValueError("time_indices must be an integer, list of integers, or None")
    
    time_indices = np.array(time_indices, dtype=int)
    n_times = len(time_indices)
    
    # Min y max global para escala de colores común
    global_vmin = min([X[species][time_indices].min() for species in species_names])
    global_vmax = max([X[species][time_indices].max() for species in species_names])
    
    # Crear figura con espacio para la barra de color global
    fig = plt.figure(figsize=(figsize_multiplier * n_times + 1, figsize_multiplier * n_species))
    
    # Crear grid con espacio para la barra de color global a la derecha
    width_ratios = [1] * n_times
    gs = fig.add_gridspec(n_species, n_times, width_ratios=width_ratios, 
                         hspace=0.4, wspace=0.3, right=0.85)  # Dejar espacio a la derecha para la barra
    
    # Lista para almacenar todas las imágenes (necesario para la barra de color)
    all_imgs = []
    
    for i, species in enumerate(species_names):
        data = X[species]
        
        for j, time_idx in enumerate(time_indices):
            ax = fig.add_subplot(gs[i, j])
            heat_data = data[time_idx]
            
            img = ax.imshow(heat_data, cmap=cmap, vmin=global_vmin, vmax=global_vmax, origin='upper', aspect='auto')
            all_imgs.append(img)  # Guardar referencia para la barra de color
            
            # Título de tiempos (solo en la primera fila)
            if i == 0:
                ax.set_title(f"t = {t[time_idx]:.1f}", fontweight='bold', pad=15)
            
            # Nombre de la especie al inicio de la fila
            if j == 0:
                ax.annotate(species, xy=(-0.1, 0.5), xycoords='axes fraction',
                            ha='right', va='center', fontweight='bold', rotation=90)
            
            # Ejes
            if i == n_species - 1:
                ax.set_xlabel(" ")
            
            ax.set_ylabel(" ")
            
            # Grilla y ticks
            rows, cols = heat_data.shape
            ax.set_xticks(np.arange(cols))
            ax.set_yticks(np.arange(rows))
            ax.grid(False)
            
            # Interacción con cursor
            cursor = mplcursors.cursor(img, hover=False)
            
            def make_callback(data_matrix):
                def on_add(sel):
                    x, y = int(sel.target[0]), int(sel.target[1])
                    if 0 <= y < data_matrix.shape[0] and 0 <= x < data_matrix.shape[1]:
                        sel.annotation.set_text(f"{data_matrix[y, x]:.3f}")
                return on_add
            
            cursor.connect("add", make_callback(heat_data))
    
    # Añadir barra de color global
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    fig.colorbar(all_imgs[0], cax=cbar_ax, label=bar_label)
    
    # Título general
    plt.suptitle(main_title, fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92, right=0.85)  # Ajustar para dejar espacio para la barra
    
    # Save the figure if requested
    if save_figure:
        # Create the directory if it doesn't exist
        os.makedirs("visualizations/plot_heatmaps_all_species_2D", exist_ok=True)
        
        # Construct the full file path
        filepath = os.path.join("visualizations", "plot_heatmaps_all_species_2D", filename)
        
        # Save the figure with high quality
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filepath}")
    
    plt.show()

# def plot_heatmaps_for_species_2D(X, species_name, t, time_indices=None, title="Heatmap for Specie"):
#     data = np.array(X[species_name])  # shape: (time, rows, cols)
#     n_time, rows, cols = data.shape
#     if time_indices is None:
#         time_indices = np.linspace(0, len(t)-1, 3, dtype=int)

#     n = len(time_indices)
#     fig = plt.figure(figsize=(5 * n, 4))
#     gs = fig.add_gridspec(1, n + 1, width_ratios=[1]*n + [0.05], wspace=0.3)
    
#     vmin = np.min([X[species_name][idx].min() for idx in time_indices])
#     vmax = np.max([X[species_name][idx].max() for idx in time_indices])
    
#     axes = [fig.add_subplot(gs[0, i]) for i in range(n)]
#     cbar_ax = fig.add_subplot(gs[0, -1])
#     imgs = []

#     for ax, idx in zip(axes, time_indices):
#         heat_data = X[species_name][idx]
#         img = ax.imshow(heat_data, cmap='viridis', vmin=vmin, vmax=vmax, origin='upper')
#         imgs.append(img)
#         ax.set_title(f"t = {t[idx]:.2f}")
#         ax.set_xlabel("") # ("Column")
#         ax.set_ylabel("") # ("Row")
#         ax.set_xticks(np.arange(cols))
#         ax.set_yticks(np.arange(rows))

#         # Cursor y función propia para cada imagen
#         cursor = mplcursors.cursor(img, hover=False)

#         def make_callback(data):
#             def on_add(sel):
#                 x, y = int(sel.target[0]), int(sel.target[1])
#                 if 0 <= y < data.shape[0] and 0 <= x < data.shape[1]:
#                     sel.annotation.set_text(f"{data[y, x]:.3f}")
#             return on_add

#         cursor.connect("add", make_callback(heat_data))

#     # Barra de color común
#     fig.colorbar(imgs[0], cax=cbar_ax)
#     plt.suptitle(f"{title} {species_name}", fontsize=16)
#     plt.tight_layout()
#     plt.show()

# Función para animar los mapas de calor con controles interactivos en html como video 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# def animate_diffusion_heatmaps_all_species_2D(t, X, species_names=None, main_title="Heatmaps for All Species", slider_label="Time", bar_label="Concentration"):
#     """
#     Animate the evolution of multiple species in heatmaps arranged by columns and save as HTML.
    
#     Parameters
#     ----------
#     t : array-like
#         Time vector    
#     X : dict
#         Dictionary {species_name: 3D array (time, rows, cols)}
#     species_names : list, optional
#         List of species to plot. If None, all species in X are plotted.
#     main_title : str
#         Main title of the figure. Default is "Heatmaps for All Species".
#     slider_label : str
#         Label for the time slider. Default is "Time".
#     bar_label : str
#         Label for the color bar. Default is "Concentration".
#     """
#     t = np.round(t, 2)

#     if species_names is None:
#         species_names = list(X.keys())
#     n_species = len(species_names)
    
#     # Calcular número de filas y columnas (máximo 5 columnas)
#     max_cols = 5
#     n_cols = min(n_species, max_cols)
#     n_rows = (n_species + n_cols - 1) // n_cols  # División entera hacia arriba

#     # Min/max global para todas las especies
#     vmin = min(np.min(X[sp]) for sp in species_names)
#     vmax = max(np.max(X[sp]) for sp in species_names)

#     # Figura y ejes (disposición por columnas)
#     fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    
#     # Si hay solo una fila o columna, convertir axes a array 2D para consistencia
#     if n_rows == 1 and n_cols == 1:
#         axes = np.array([[axes]])
#     elif n_rows == 1:
#         axes = axes.reshape(1, -1)
#     elif n_cols == 1:
#         axes = axes.reshape(-1, 1)
    
#     # Ajustar espacio
#     plt.subplots_adjust(bottom=0.1, top=0.9, hspace=0.4, wspace=0.3)

#     ims = []
#     for i, sp in enumerate(species_names):
#         row_idx = i // n_cols
#         col_idx = i % n_cols
#         ax = axes[row_idx, col_idx]
        
#         data = np.array(X[sp])
#         im = ax.imshow(data[0], cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
#         ax.set_title(f"{sp} (t = {t[0]:.2f})", fontweight="bold")
#         ax.set_xlabel(" ") # Column")
#         ax.set_ylabel(" ") # Row")
#         rows, cols = data.shape[1:]
#         ax.set_xticks(np.arange(cols))
#         ax.set_yticks(np.arange(rows))
#         ims.append((im, data, ax, sp))

#     # Ocultar ejes vacíos si el número de especies no llena completamente la cuadrícula
#     total_cells = n_rows * n_cols
#     if total_cells > n_species:
#         for i in range(n_species, total_cells):
#             row_idx = i // n_cols
#             col_idx = i % n_cols
#             axes[row_idx, col_idx].axis('off')

#     # Un único colorbar global
#     cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # [left, bottom, width, height]
#     fig.colorbar(ims[0][0], cax=cbar_ax, label=bar_label)

#     def update(frame):
#         for im, data, ax, sp in ims:
#             im.set_data(data[frame])
#             ax.set_title(f"{sp} (t = {t[frame]:.2f})", fontweight="bold")
#         return [item[0] for item in ims]

#     ani = FuncAnimation(fig, update, frames=len(t), interval=50, blit=False)

#     plt.suptitle(main_title, fontsize=16, fontweight="bold", y=0.97)

#     # Generar HTML
#     html = ani.to_html5_video()

#     # Crear carpeta si no existe
#     folder_path = 'visualizations/plot_heatmaps_all_species_2D_animate'
#     os.makedirs(folder_path, exist_ok=True)

#     # Guardar en archivo HTML
#     file_path = os.path.join(folder_path, 'animation.html')
#     with open(file_path, 'w') as f:
#         f.write(html)
#     print(f"Animation saved as: {file_path}")


# Función para animar los mapas de calor con controles interactivos
def animate_diffusion_heatmaps_all_species_2D(t, X, species_names=None, main_title="Heatmaps for All Species", slider_label="Time", bar_label="Concentration"):
    """
    Animate the evolution of multiple species in heatmaps arranged by columns.
    
    Parameters
    ----------
    t : array-like
        Time vector    
    X : dict
        Dictionary {species_name: 3D array (time, rows, cols)}
    species_names : list, optional
        List of species to plot. If None, all species in X are plotted.
    main_title : str
        Main title of the figure. Default is "Heatmaps for All Species".
    slider_label : str
        Label for the time slider. Default is "Time".
    bar_label : str
        Label for the color bar. Default is "Concentration".
    """
    t = np.round(t, 2)

    if species_names is None:
        species_names = list(X.keys())
    n_species = len(species_names)
    
    # Calcular número de filas y columnas (máximo 5 columnas)
    max_cols = 5
    n_cols = min(n_species, max_cols)
    n_rows = (n_species + n_cols - 1) // n_cols  # División entera hacia arriba

    # Min/max global para todas las especies
    vmin = min(np.min(X[sp]) for sp in species_names)
    vmax = max(np.max(X[sp]) for sp in species_names)

    # Figura y ejes (disposición por columnas)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    
    # Si hay solo una fila o columna, convertir axes a array 2D para consistencia
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Ajustar espacio para el slider y botones
    plt.subplots_adjust(bottom=0.25, top=0.9, hspace=0.4, wspace=0.3)

    ims = []
    for i, sp in enumerate(species_names):
        row_idx = i // n_cols
        col_idx = i % n_cols
        ax = axes[row_idx, col_idx]
        
        data = np.array(X[sp])
        im = ax.imshow(data[0], cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
        ax.set_title(f"{sp}", fontweight="bold") # (t = {t[0]:.2f})", fontweight="bold")
        ax.set_xlabel(" ") # Column")
        ax.set_ylabel(" ") # Row")
        rows, cols = data.shape[1:]
        ax.set_xticks(np.arange(cols))
        ax.set_yticks(np.arange(rows))
        ims.append((im, data, ax))

    # Ocultar ejes vacíos si el número de especies no llena completamente la cuadrícula
    total_cells = n_rows * n_cols
    if total_cells > n_species:
        for i in range(n_species, total_cells):
            row_idx = i // n_cols
            col_idx = i % n_cols
            axes[row_idx, col_idx].axis('off')

    # Un único colorbar global
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.6])  # [left, bottom, width, height]
    fig.colorbar(ims[0][0], cax=cbar_ax, label=bar_label)

    # Slider
    ax_slider = plt.axes([0.2, 0.12, 0.6, 0.03])
    slider = Slider(ax_slider, slider_label, t[0], t[-1], valinit=t[0], valstep=(t[1]-t[0]))

    # Botones
    ax_play = plt.axes([0.3, 0.05, 0.1, 0.05])
    ax_pause = plt.axes([0.6, 0.05, 0.1, 0.05])
    btn_play = Button(ax_play, 'Play')
    btn_pause = Button(ax_pause, 'Pause')

    playing = [False]

    def update(val):
        idx = np.argmin(np.abs(t - val))
        for im, data, ax, sp in zip([item[0] for item in ims], 
                                   [item[1] for item in ims], 
                                   [item[2] for item in ims], 
                                   species_names):
            im.set_data(data[idx])
            ax.set_title(f"{sp}", fontweight="bold") # (t = {t[idx]:.2f})", fontweight="bold")
        fig.canvas.draw_idle()

    def play(event):
        playing[0] = True

    def pause(event):
        playing[0] = False

    btn_play.on_clicked(play)
    btn_pause.on_clicked(pause)
    slider.on_changed(update)

    def animate(frame):
        if playing[0]:
            current_val = slider.val
            next_val = current_val + (t[1] - t[0])
            if next_val > t[-1]:
                playing[0] = False  # Detener la animación al llegar al final
            else:
                slider.set_val(next_val)

    ani = animation.FuncAnimation(fig, animate, interval=50)
    update(t[0])  # inicializa el primer frame

    plt.suptitle(main_title, fontsize=16, fontweight="bold", y=0.97)
    plt.show()

########################################################################################

def plot_species_dynamics_MP(t, time_series, species, num_patches=4, separate_plots=False, 
                         filename='dynamics.png', title='Dinámicas de Especies en Cada Parche', 
                         figsize=(12, 8)):
    """
    Grafica las series temporales de las concentraciones de especies en cada parche.
    
    Parámetros:
    - t: Arreglo de tiempos (numpy array).
    - time_series: DataFrame con las concentraciones, columnas en formato (especie, parche).
    - species: Lista de nombres de las especies (ej. ['l', 's1', 's2']).
    - num_patches: Número de parches (default: 4).
    - separate_plots: Si True, crea un subgráfico por especie; si False, grafica todo en una figura (default: False).
    - filename: Nombre del archivo para guardar el gráfico (default: 'dynamics.png').
    - title: Título del gráfico (default: 'Dinámicas de Especies en Cada Parche').
    - figsize: Tupla con el tamaño de la figura (default: (12, 8)).
    """
    if separate_plots:
        # Crear subgráficos para cada especie
        fig, axs = plt.subplots(len(species), 1, figsize=(figsize[0], figsize[1] * len(species) / 2), sharex=True)
        for i, sp in enumerate(species):
            for p in range(num_patches):
                axs[i].plot(t, time_series[(sp, p)], label=f'Parche {p}')
            axs[i].set_title(f'Concentración de {sp}')
            axs[i].set_ylabel('Concentración')
            axs[i].legend()
            axs[i].grid(True)
        axs[-1].set_xlabel('Tiempo')
        plt.suptitle(title, y=1.02)
    else:
        # Graficar todas las especies en una sola figura
        plt.figure(figsize=figsize)
        for sp in species:
            for p in range(num_patches):
                plt.plot(t, time_series[(sp, p)], label=f'{sp}, Parche {p}')
        plt.xlabel('Tiempo')
        plt.ylabel('Concentración')
        plt.title(title)
        plt.legend()
        plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

import matplotlib.pyplot as plt
import numpy as np
import mplcursors

def plot_heatmaps_for_species_2d(time_series, species, t, patch_shape=(2, 2), time_indices=None, 
                                filename_prefix='heatmap_2d', figsize=(5, 4)):
    """
    Grafica mapas de calor 2D para cada especie en una cuadrícula de parches en múltiples puntos temporales.
    
    Parámetros:
    - time_series: DataFrame con las concentraciones, columnas en formato (especie, parche).
    - species: Lista de nombres de las especies (ej. ['l', 's1', 's2']).
    - t: Arreglo de tiempos (numpy array).
    - patch_shape: Tupla (filas, columnas) de la cuadrícula de parches (default: (2, 2)).
    - time_indices: Lista de índices de tiempo a graficar (default: None, usa 3 tiempos equiespaciados).
    - filename_prefix: Prefijo para los nombres de los archivos guardados (default: 'heatmap_2d').
    - figsize: Tupla con el tamaño base de cada subgráfico (default: (5, 4)).
    """
    rows, cols = patch_shape
    num_patches = rows * cols
    num_species = len(species)
    
    if time_indices is None:
        time_indices = np.linspace(0, len(t)-1, 3, dtype=int)  # Inicio, medio, final
    
    n_times = len(time_indices)
    
    # Crear una figura por especie
    for sp in species:
        # Obtener los valores mínimo y máximo para la especie en los tiempos seleccionados
        vmin = np.min([time_series[(sp, p)].iloc[time_indices] for p in range(num_patches)])
        vmax = np.max([time_series[(sp, p)].iloc[time_indices] for p in range(num_patches)])
        
        # Configurar la figura
        fig = plt.figure(figsize=(figsize[0] * n_times, figsize[1]))
        gs = fig.add_gridspec(1, n_times + 1, width_ratios=[1]*n_times + [0.05], wspace=0.3)
        axes = [fig.add_subplot(gs[0, i]) for i in range(n_times)]
        cbar_ax = fig.add_subplot(gs[0, -1])
        imgs = []
        
        # Graficar un mapa de calor para cada tiempo seleccionado
        for ax, idx in zip(axes, time_indices):
            # Obtener concentraciones para el tiempo idx
            conc_data = np.array([time_series[(sp, p)].iloc[idx] for p in range(num_patches)]).reshape(rows, cols)
            
            # Mapa de calor
            img = ax.imshow(conc_data, cmap='viridis', vmin=vmin, vmax=vmax, origin='upper')
            imgs.append(img)
            ax.set_title(f"t = {t[idx]:.2f}")
            ax.set_xlabel("Columna")
            ax.set_ylabel("Fila")
            ax.set_xticks(np.arange(cols))
            ax.set_yticks(np.arange(rows))
            
            # Añadir valores numéricos
            for r in range(rows):
                for c in range(cols):
                    ax.text(c, r, f'{conc_data[r, c]:.2f}', ha='center', va='center', color='white')
            
            # Añadir interacción con mplcursors
            cursor = mplcursors.cursor(img, hover=False)
            
            def make_callback(data):
                def on_add(sel):
                    x, y = int(sel.target[0]), int(sel.target[1])
                    if 0 <= y < data.shape[0] and 0 <= x < data.shape[1]:
                        sel.annotation.set_text(f"{data[y, x]:.3f}")
                return on_add
            
            cursor.connect("add", make_callback(conc_data))
        
        # Añadir barra de color común
        fig.colorbar(imgs[0], cax=cbar_ax)
        plt.suptitle(f"Especie {sp}", fontsize=16)
        plt.tight_layout()
        plt.savefig(f"{filename_prefix}_{sp}.png")
        plt.show()

# Función para animar mapas de calor 2D
def animate_diffusion_heatmaps_for_species_2d(time_series, species_name, t, patch_shape=(2, 2), 
                                             filename=None, figsize=(6, 5)):
    """
    Anima mapas de calor 2D para una especie en una cuadrícula de parches con controles interactivos.
    """
    rows, cols = patch_shape
    num_patches = rows * cols
    
    # Redondear tiempos a dos decimales
    t = np.round(t, 2)
    
    # Extraer datos para la especie específica y reorganizar en (time, rows, cols)
    data = np.array([time_series[(species_name, p)] for p in range(num_patches)]).T.reshape(-1, rows, cols)
    n_time, _, _ = data.shape
    
    # Calcular valores mínimo y máximo para la escala de color
    vmin = np.min(data)
    vmax = np.max(data)
    
    # Configurar la figura
    fig, ax = plt.subplots(figsize=figsize)
    plt.subplots_adjust(bottom=0.35)  # Espacio para controles
    
    # Mapa de calor inicial
    im = ax.imshow(data[0], cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
    ax.set_title(f"Especie {species_name} (t = {t[0]:.2f})")
    ax.set_xlabel("Columna")
    ax.set_ylabel("Fila")
    ax.set_xticks(np.arange(cols))
    ax.set_yticks(np.arange(rows))
    
    # Añadir valores numéricos en cada celda
    texts = []
    for r in range(rows):
        for c in range(cols):
            text = ax.text(c, r, f'{data[0, r, c]:.2f}', ha='center', va='center', color='white')
            texts.append(text)
    
    # Barra de color
    cbar = fig.colorbar(im, ax=ax, orientation='vertical')
    
    # Slider
    ax_slider = plt.axes([0.2, 0.2, 0.6, 0.03])
    slider = Slider(ax_slider, 'Time', t[0], t[-1], valinit=t[0], valstep=(t[1] - t[0]))
    
    # Botones
    ax_play = plt.axes([0.3, 0.1, 0.1, 0.05])
    ax_pause = plt.axes([0.6, 0.1, 0.1, 0.05])
    btn_play = Button(ax_play, 'Play')
    btn_pause = Button(ax_pause, 'Pause')
    
    playing = [False]  # Estado mutable para la animación
    
    def update(val):
        idx = np.argmin(np.abs(t - val))
        im.set_data(data[idx])
        ax.set_title(f"Especie {species_name} (t = {t[idx]:.2f})")
        # Actualizar valores numéricos
        for r in range(rows):
            for c in range(cols):
                texts[r * cols + c].set_text(f'{data[idx, r, c]:.2f}')
        fig.canvas.draw_idle()
    
    def play(event):
        playing[0] = True
    
    def pause(event):
        playing[0] = False
    
    btn_play.on_clicked(play)
    btn_pause.on_clicked(pause)
    slider.on_changed(update)
    
    def animate(frame):
        if playing[0]:
            current_val = slider.val
            next_val = current_val + (t[1] - t[0])
            if next_val > t[-1]:
                next_val = t[0]  # Loop
            slider.set_val(next_val)
    
    ani = animation.FuncAnimation(fig, animate, interval=20, frames=n_time)
    
    update(t[0])  # Inicializa el primer frame
    
    if filename:
        # Guardar la animación como MP4 (requiere ffmpeg)
        ani.save(filename, writer='ffmpeg', fps=50)
    
    plt.show()
################################################################################
# Plot the number of abstractions over time
################################################################################
import matplotlib.pyplot as plt
import os

def plot_abstraction_size(abstract_time_series, xlabel="Time", ylabel="Number of Species", 
                         title="Number of species per abstraction over time", marker='o', 
                         label="Abstraction Size", save_figure=True, filename="abstraction_size_plot.png"):
    """
    Plots the number of abstractions over the time series and optionally saves the figure.

    Parameters:
    abstract_time_series (pd.DataFrame): Time series with a 'Time' column and a column for species abstractions.
    xlabel (str): Label for the x-axis. Default is "Time".
    ylabel (str): Label for the y-axis. Default is "Number of Species".
    title (str): Title of the plot. Default is "Number of species per abstraction over time".
    marker (str): Marker style for the plot. Default is 'o'.
    label (str): Legend label for the plot. Default is "Abstraction Size".
    save_figure (bool): Whether to save the figure. Default is True.
    filename (str): Name of the file to save. Default is "abstraction_size_plot.png".

    Raises:
    ValueError: If the DataFrame does not contain a 'Time' column.

    Returns:
    None: Displays a line plot of abstraction sizes over time.
    """
    if 'Time' not in abstract_time_series.columns:
        raise ValueError("The DataFrame must include a 'Time' column for time values.")
    
    # Create a new figure and axis object for the plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Plot size 10x6 inches
    
    # Create a mapping of time and number of active species (abstraction size)
    abstraction_mapping = {
        idx: (row["Time"], len(row["Abstraction"]))  # Map index to (time, abstraction size)
        for idx, row in abstract_time_series.iterrows()
    }
    
    # Extract times and abstraction sizes
    times = [abstraction_mapping[idx][0] for idx in abstraction_mapping]  # Extract time values
    sizes = [abstraction_mapping[idx][1] for idx in abstraction_mapping]  # Extract abstraction sizes

    # Plot abstractions as ordered points
    ax.plot(times, sizes, marker=marker, linestyle='-', label=label)  # Plot with line and marker
    ax.set_xlabel(xlabel)  # Set x-axis label
    ax.set_ylabel(ylabel)  # Set y-axis label
    ax.set_title(title)  # Set plot title
    ax.grid()  # Enable grid on the plot
    ax.legend()  # Add a legend for the abstraction size
    
    plt.tight_layout()  # Adjust layout to fit all elements
    
    # Save the figure if requested
    if save_figure:
        # Create the directory if it doesn't exist
        os.makedirs("visualizations/plot_abstraction_size", exist_ok=True)
        
        # Construct the full file path
        filepath = os.path.join("visualizations", "plot_abstraction_size", filename)
        
        # Save the figure with high quality
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filepath}")
    
    plt.show()  # Display the plot

######################################################################################## 

import matplotlib.pyplot as plt
import pandas as pd
import os

def plot_abstraction_sets(abstract_time_series, xlabel="Time", ylabel="Species", 
                         title="Abstraction of Time Series", save_figure=True, 
                         filename="abstraction_sets_plot.png"):
    """
    Plots the time series abstraction generated by the abstraction_ordinary function and optionally saves the figure.

    Parameters:
    abstract_time_series (pd.DataFrame): DataFrame containing 'Time' and 'Abstraction' columns.
                                         'Abstraction' should be a list of species present at each time point.
    xlabel (str): Label for the x-axis. Default is "Time".
    ylabel (str): Label for the y-axis. Default is "Presence of Species".
    title (str): Title of the plot. Default is "Abstraction of ODE Time Series".
    save_figure (bool): Whether to save the figure. Default is True.
    filename (str): Name of the file to save. Default is "abstraction_sets_plot.png".

    Returns:
    None: Displays a stacked bar chart showing the presence of species over time.
    """
    # Extract all unique species from the 'Abstraction' column
    all_species = sorted({species for row in abstract_time_series["Abstraction"] for species in row})  # Unique sorted species
    
    # Build a binary matrix indicating the presence of species over time
    binary_matrix = pd.DataFrame(
        [
            {species: (species in row) for species in all_species}  # Check presence of each species in a row
            for row in abstract_time_series["Abstraction"]
        ],
        index=abstract_time_series["Time"]  # Use time as the index
    ).astype(int)  # Convert to integer (0 or 1)
    
    # Create a new figure and axis object for the plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Plot size 10x6 inches
    
    # Plot each species as a stacked bar chart
    for i, species in enumerate(all_species):
        ax.bar(
            binary_matrix.index,  # Time points
            binary_matrix[species],  # Presence of species
            label=species,  # Label for legend
            bottom=binary_matrix.iloc[:, :i].sum(axis=1)  # Stacked position
        )
    
    # Configure plot labels and title
    ax.set_xlabel(xlabel)  # Set x-axis label
    ax.set_ylabel(ylabel)  # Set y-axis label
    ax.set_title(title)  # Set plot title
    ax.legend(title="Species")  # Add legend with title
    
    # Save the figure if requested
    if save_figure:
        # Create the directory if it doesn't exist
        os.makedirs("visualizations/plot_abstraction_sets", exist_ok=True)
        
        # Construct the full file path
        filepath = os.path.join("visualizations", "plot_abstraction_sets", filename)
        
        # Save the figure with high quality
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filepath}")
    
    plt.show()  # Display the plot

######################################################################################
# Abstraction graph static and movie 
######################################################################################
# # # Original: plot_static_abstraction_graph
def plot_static_abstraction_graph(abstract_time_series, title="Static Abstraction Graph"):
    """
    Plots a static abstraction graph with nodes, edges, and associated attributes.

    The graph represents:
    - Node size: Depends on the frequency of occurrences of each abstraction.
    - Node color: Represents the weighted average time of occurrence (normalized).
    - Edge thickness: Indicates the frequency of transitions between abstractions.

    Parameters:
    abstract_time_series (pd.DataFrame): A DataFrame with the following columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    title (str): Title of the graph. Defaults to "Static Abstraction Graph".

    Returns:
    None. Displays a static plot of the abstraction graph.
    """
    abstractions = abstract_time_series["Abstraction"]  # Extracts the abstraction column
    times = abstract_time_series["Time"]  # Extracts the time column
    nodes = abstractions.apply(tuple).value_counts()  # Counts the occurrences of each abstraction

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  # Transitions between consecutive abstractions
    transitions_freq = pd.Series(transitions).value_counts()  # Counts the frequency of each transition

    # Create the graph
    G = nx.DiGraph()  # Creates a directed graph

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]  # Filters the times associated with the current node
        weighted_time = (node_times * freq).mean()  # Calculates the weighted average time
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())  # Normalizes the time
        # color = plt.cm.coolwarm(1 - normalized_time)  # Assigns a color based on the normalized time
        color = to_hex(plt.cm.coolwarm(1 - normalized_time))  # Convierte el color a formato hexadecimal
        G.add_node(node, size=freq, color=color)  # Adds the node to the graph with attributes

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)  # Adds edges to the graph with the corresponding weight

    # Hierarchical layout
    pos = nx.planar_layout(G) # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")  # Calculates the node positions in a hierarchical layout
    pos = {node: (x, -y) for node, (x, y) in pos.items()}  # Inverts the y-axis so the graph is drawn from bottom to top

    # Create the plot
    plt.figure(figsize=(12, 8))  # Creates a figure with a specific size

    # Draw nodes
    node_sizes = [G.nodes[node]["size"] * 100 for node in G.nodes]  # Scales the node sizes
    node_colors = [G.nodes[node]["color"] for node in G.nodes]  # Gets the node colors
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors)  # Draws the nodes

    # Draw edges
    edge_weights = [G.edges[edge]["weight"] for edge in G.edges]  # Gets the edge weights
    nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.7)  # Draws the edges with adjusted opacity

    # Draw node labels
    # nx.draw_networkx_labels(G, pos, font_size=8, font_color="black")  # Draws node labels 
    node_labels = {node: f"({', '.join(map(str, node))})" for node in G.nodes}  # Formatea las etiquetas
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_color="black")  # Usa las etiquetas formateadas

    # Title and display
    plt.title(title)  # Adds a title to the graph
    plt.axis("off")  # Hides the axes
    plt.show()  # Displays the graph

########################################################################################

def get_plot_static_abstraction_graph_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph.html"
):
    """
    Generates an interactive HTML visualization of a directed graph.

    The graph is constructed from a time series of abstractions, with nodes representing unique abstractions
    and edges representing transitions between consecutive abstractions. Node attributes are customized
    based on their frequency, and edge attributes are customized based on transition frequencies.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Abstraction": A column with abstraction data (e.g., lists or tuples).
        - "Time": A column with corresponding time values for each abstraction.
    node_size : int, optional
        The size of the nodes in the visualization (default is 20).
    edge_color : str, optional
        The color of the edges in the graph (default is "gray").
    shape_node : str, optional
        The shape of the nodes (e.g., 'dot', 'square'; default is 'dot').
    node_font_size : int, optional
        The font size for node labels (default is 14).
    edge_width : int, optional
        The base width of the edges (default is 2).
    filename : str, optional
        The name of the output HTML file (default is "Static_Abstraction_Graph.html").

    Returns:
    -------
    pyvis.network.Network
        A PyVis Network object containing the interactive visualization.

    Notes:
    -----
    - The node color encodes its position in the sequence:
      * Blue: Initial node
      * Red: Final node
      * Gradient colors (coolwarm colormap): Intermediate nodes based on normalized time.
    - Nodes and edges are added to the visualization with specific attributes such as size, color, 
      label, and weight.
    - The output file is an interactive HTML visualization saved at the specified filename.

    Example:
    -------
    >>> df = pd.DataFrame({
    >>>     "Time": [0, 1, 2, 3, 4],
    >>>     "Abstraction": [[A], [A,B], [C], [A,C], [A]]    
    >>> })
    >>> get_plot_static_abstraction_graph_html(df) 
    """
    # Extract abstraction and time columns
    abstractions = abstract_time_series["Abstraction"]  
    times = abstract_time_series["Time"]  
    
    # Count node frequencies
    nodes = abstractions.apply(tuple).value_counts(sort=False)  
    # print("Nodes de", nodes)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  
    transitions_freq = pd.Series(transitions).value_counts()  
    # print("transitions", transitions)
    # print("transitions_freq", transitions_freq)
     
    # Create the graph
    G = nx.DiGraph()  

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}
        node_times = times[abstractions.apply(tuple) == node]  
        weighted_time = (node_times * freq).mean()  
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())  

        if node == tuple(abstractions.iloc[0]):  
            color = "#0000FF"  
        elif node == tuple(abstractions.iloc[-1]):  
            color = "#FF0000"  
        else:
            color = to_hex(plt.cm.coolwarm(1 - normalized_time))  

        label = f"{node_labels[node]} ({freq})"  
        G.add_node(node, size=freq, color=color, label=label)  

    # Add edges with attributes
    for (source, target), freq in transitions_freq.items():
        G.add_edge(
            source, 
            target, 
            weight=freq,  
            color=edge_color,  
            width=edge_width + freq / max(transitions_freq) * edge_width  
        )

    # Convert to PyVis visualization 
    net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "nodes": {{
        "shape": "{shape_node}",  
        "font": {{
        "size": {node_font_size},  
        "align": "center"  
        }},
        "borderWidth": 2,  
        "borderColor": "black"  
      }},
      "edges": {{
        "smooth": false,
        "color": "{edge_color}",
        "width": {edge_width}  
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{
          "enabled": false
        }},
        "hierarchicalRepulsion": {{
          "nodeDistance": 150
        }}
      }},
      "layout": {{
        "hierarchical": {{
          "enabled": false,
          "direction": "DU",  
          "sortMethod": "directed"
        }}
      }}
    }}
    """) 

    # Add nodes to the PyVis network
    for node, attrs in G.nodes(data=True):
        net.add_node(
            node_labels[node],  
            label=attrs['label'],  
            title=f"{node}",
            color=attrs['color'],
            size=node_size
        )

    # Add edges to the PyVis network
    for source, target, attrs in G.edges(data=True):
        net.add_edge(
            node_labels[source],  
            node_labels[target],  
            width=edge_width
        )
        
    # Save the visualization to an HTML file 
    net.html = net.generate_html() 
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return net

########################################################################################

def plot_static_abstraction_graph_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph.html"
):
    """
    Generate and visualize the static abstraction graph as an interactive HTML file.

    This function is a wrapper to visualize the containment hierarchy among sets
    using a static abstraction graph. The graph is saved as an HTML file and can be 
    opened in a web browser.

    Args:
        abstract_time_series (pd.DataFrame): DataFrame containing abstraction and time series information.
            It must include at least the following columns:
            - "Abstraction": Represents the sets or elements to visualize.
            - "Time": Represents the corresponding time for the abstractions.
        node_size (int, optional): Size of the nodes in the graph. Default is 20.
        edge_color (str, optional): Color of the edges connecting nodes. Default is "gray".
        shape_node (str, optional): Shape of the nodes. Options include 'dot', 'square', 
            'triangle', 'star', or 'diamond'. Default is 'dot'.
        node_font_size (int, optional): Font size for node labels. Default is 14.
        edge_width (int, optional): Thickness of the edges. Default is 2.
        filename (str, optional): Name of the output HTML file. Default is "Static_Abstraction_Graph.html".

    Raises:
        FileNotFoundError: If the generated HTML file is not found in the specified location.

    Side Effects:
        - Saves the generated graph as an HTML file at the specified location.
        - Opens the HTML file in the default web browser.

    Example:
        ```python
        plot_static_abstraction_graph_html(
            abstract_time_series=my_dataframe,
            node_size=30,
            edge_color="blue",
            shape_node="triangle",
            filename="My_Abstraction_Graph.html"
        )
        ```
    """
    # Call the `get_plot_static_abstraction_graph_html` function to generate the graph
    get_plot_static_abstraction_graph_html(
        abstract_time_series,  
        node_size=node_size,  
        edge_color=edge_color, 
        shape_node=shape_node, 
        node_font_size=node_font_size, 
        edge_width=edge_width,  
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `get_plot_static_abstraction_graph_html` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The hierarchy visualization was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

########################################################################################

def get_plot_static_abstraction_graph_hierarchy_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph_Hierarchy.html"
):
    """
    Generates a hierarchical abstraction graph from time series data and saves it as an interactive HTML file.

    Args:
        abstract_time_series (pd.DataFrame): A DataFrame with two columns:
            - "Abstraction": Contains sets or lists representing the abstract states.
            - "Time": Contains the corresponding timestamps for each state.
        node_size (int): Size of the nodes in the graph visualization.
        edge_color (str): Color of the edges in the graph.
        shape_node (str): Shape of the nodes; options include 'dot', 'square', 'triangle', etc.
        node_font_size (int): Font size for the node labels.
        edge_width (int): Base width of the edges in the graph.
        filename (str): Name of the output HTML file where the visualization will be saved.

    Returns:
        net (pyvis.network.Network): A PyVis Network object representing the generated graph.

    Description:
        - Extracts nodes and their frequencies from the "Abstraction" column of the time series data.
        - Computes containment relationships among nodes to build a directed graph.
        - Nodes are colored based on their temporal positions (initial, final, or normalized weighted time).
        - Saves the graph as an interactive HTML file with hierarchical layout and opens it in the browser.

    Example:
        ```python
        abstract_time_series = pd.DataFrame({
            "Abstraction": [{"A", "B"}, {"A"}, {"A", "C"}, {"C"}],
            "Time": [1, 2, 3, 4]
        })
        get_plot_static_abstraction_graph_hierarchy_html(
            abstract_time_series,
            node_size=30,
            edge_color="blue",
            shape_node="square",
            filename="Hierarchy_Graph.html"
        )
        ```
    """
    abstractions = abstract_time_series["Abstraction"]  # Extracts the abstraction column
    times = abstract_time_series["Time"]  # Extracts the time column
    nodes = abstractions.apply(tuple).value_counts(sort=False)  # Counts the occurrences of each abstraction    
    # print("nodes", nodes)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  # Transitions between consecutive abstractions
    transitions_freq = pd.Series(transitions).value_counts()  # Counts the frequency of each transition
    # print("transitions", transitions)
    # print("transitions_freq", transitions_freq)
    
    # Build the containment graph
    containment_graph = {node: [] for node in nodes.index}  # Initialize containment graph
    for i, child_set in enumerate(nodes.index):  # Iterate over sets as children
        for j, parent_set in enumerate(nodes.index):  # Iterate over sets as parents
            if i != j and set(child_set).issubset(parent_set):  # If they are not the same and containment exists
                if parent_set not in containment_graph[child_set]:  # Check subsets of child_set
                    containment_graph[child_set].append(parent_set)  # Add the relationship
    # print("containment_graph", containment_graph)

    # Create the graph
    G = nx.DiGraph()  # Directed graph

    # Create unique labels for nodes (e.g., S1, S2, S3...)
    node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}

    # Add nodes with attributes 
    for node, freq in nodes.items():  
        node_times = times[abstractions.apply(tuple) == node]  # Filters the times associated with the current node
        weighted_time = (node_times * freq).mean()  # Calculates the weighted average time
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())  # Normalizes the time

        if node == tuple(abstractions.iloc[0]):  # Initial node
            color = "#0000FF"  # Blue
        elif node == tuple(abstractions.iloc[-1]):  # Final node
            color = "#FF0000"  # Red
        else:
            color = to_hex(plt.cm.coolwarm(1 - normalized_time))  # Color based on normalized time

        label = f"{node_labels[node]} ({freq})"  # Label includes the unique label, frequency, and elements
        G.add_node(node, size=freq, color=color, label=label)  # Add node with attributes

    # print("nodes items:", nodes)

    # Add edges based on containment graph (parents -> children) 
    for child, parents in containment_graph.items(): 
        for parent in parents:
            G.add_edge(child, parent, weight=1)  # Add directed edges from parent to child

    # Initialize the PyVis network
    net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "nodes": {{
        "shape": "{shape_node}",  
        "font": {{
          "size": {node_font_size},  
          "align": "center"  
        }},
        "borderWidth": 2,  
        "borderColor": "black"  
      }},
      "edges": {{
        "smooth": false,
        "color": "{edge_color}",
        "width": {edge_width}  
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{
          "enabled": false
        }},
        "hierarchicalRepulsion": {{
          "nodeDistance": 150
        }}
      }},
      "layout": {{
        "hierarchical": {{
          "enabled": true,
          "direction": "DU",  
          "sortMethod": "directed"
        }}
      }}
    }}
    """)

    # Add nodes to the PyVis network
    for node, attrs in G.nodes(data=True):
        net.add_node(
            node_labels[node],  # Use the unique label for the node
            label=attrs['label'],  # Set the label with frequency and elements
            title=f"{node}",
            color=attrs['color'],
            size=node_size
        )

    # Add edges to the PyVis network
    for source, target, attrs in G.edges(data=True):
        net.add_edge(
            node_labels[source],  # Use unique labels for source
            node_labels[target],  # Use unique labels for target
            width=edge_width
        )

    # Save the visualization to an HTML file 
    net.html = net.generate_html() 
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return net

########################################################################################

def plot_static_abstraction_graph_hierarchy_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph_Hierarchy.html"
):        
    """
    Generates and visualizes the containment hierarchy of abstractions as an interactive HTML file.
    
    This function serves as a wrapper to call the underlying function 
    `get_plot_static_abstraction_graph_hierarchy_html`, which creates the graph visualization. 
    Once the graph is generated, it saves the visualization to an HTML file and optionally opens it 
    in the user's default web browser.
    
    Args:
        abstract_time_series (pd.DataFrame): A DataFrame containing columns "Abstraction" (list or set of 
            elements) and "Time" (time associated with the abstractions). 
        node_size (int, optional): Size of the nodes in the visualization. Default is 20.
        edge_color (str, optional): Color of the edges in the graph. Default is "gray".
        shape_node (str, optional): Shape of the nodes in the graph. Options include 'dot', 'square', 
            'triangle', 'star', and 'diamond'. Default is 'dot'.
        node_font_size (int, optional): Font size of the labels displayed on nodes. Default is 14.
        edge_width (int, optional): Width of the edges connecting nodes. Default is 2.
        filename (str, optional): Name of the output HTML file to save the graph visualization. Default is 
            "Static_Abstraction_Graph_Hierarchy.html".
    
    Raises:
        FileNotFoundError: If the HTML file is not found after attempting to generate it.
    
    Outputs:
        The function generates an interactive HTML file that represents the abstraction graph hierarchy 
        and opens it in the default web browser.
    
    Example:
        >>> plot_static_abstraction_graph_hierarchy_html(
        ...     abstract_time_series=data,
        ...     node_size=30,
        ...     edge_color="blue",
        ...     shape_node="triangle",
        ...     filename="AbstractionGraph.html"
        ... )
    """
    # Call the function to generate the visualization
    get_plot_static_abstraction_graph_hierarchy_html(
        abstract_time_series,  
        node_size=node_size,  
        edge_color=edge_color, 
        shape_node=shape_node, 
        node_font_size=node_font_size, 
        edge_width=edge_width,  
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `get_plot_static_abstraction_graph_hierarchy_html` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The hierarchy visualization was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

########################################################################################

def plot_static_abstraction_graph_hierarchy(abstract_time_series, title="Hierarchical Abstraction Graph"):
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    nodes = abstractions.apply(tuple).value_counts()

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.coolwarm(1 - normalized_time))
        G.add_node(node, size=freq, color=color)

    # Add edges for transitions
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Add edges based on containment relationships
    node_list = list(G.nodes)
    for i, smaller in enumerate(node_list):
        for j, larger in enumerate(node_list):
            if i != j and set(smaller).issubset(set(larger)):  # Containment relationship
                if not G.has_edge(smaller, larger):  # Avoid duplicate edges
                    G.add_edge(smaller, larger, weight=1)  # Default weight

    # Compute hierarchical positions based on containment
    levels = {}
    for node in G.nodes:
        level = sum(1 for other in G.nodes if set(node).issubset(set(other)) and node != other)
        levels[node] = level

    pos = {node: (0, -level) for node, level in levels.items()}

    # Create the plot
    plt.figure(figsize=(12, 8))

    # Draw nodes
    node_sizes = [G.nodes[node]["size"] * 100 for node in G.nodes]
    node_colors = [G.nodes[node]["color"] for node in G.nodes]
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors)

    # Draw edges
    edge_weights = [G.edges[edge].get("weight", 1) for edge in G.edges]  # Use default weight
    nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.7)

    # Draw node labels
    node_labels = {node: f"({', '.join(map(str, node))})" for node in G.nodes}
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_color="black")

    # Title and display
    plt.title(title)
    plt.axis("off")
    plt.show()

########################################################################################

def plot_static_abstraction_graph_html_with_hierarchy(
    abstract_time_series,
    filename="static_abstraction_graph_hierarchy.html",
    node_size_scale=5,
    edge_width_scale=0.5,
    default_color="cyan",
    title="plot_static_abstraction_graph_html_with_hierarchy"
):
    """
    Generates an HTML visualization of a static abstraction graph considering hierarchy based on subset relationships.

    Parameters:
    abstract_time_series (pd.DataFrame): DataFrame with columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    filename (str): Name of the output HTML file. Defaults to "static_abstraction_graph_hierarchy.html".
    node_size_scale (float): Scaling factor for node sizes.
    edge_width_scale (float): Scaling factor for edge widths.
    default_color (str): Default color for nodes.
    title (str): Title of the graph. Defaults to "Static Abstraction Graph with Hierarchy".

    Returns:
    str: Filename of the generated HTML file.
    """
    # Extract data from DataFrame
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    nodes = abstractions.apply(tuple).value_counts()

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create hierarchical levels based on subset relationships
    hierarchy = {}
    sorted_nodes = sorted(nodes.keys(), key=lambda x: (len(x), x))  # Sort by size of the set and lexicographically
    for node in sorted_nodes:
        hierarchy[node] = []
        for potential_parent in sorted_nodes:
            if set(node).issubset(set(potential_parent)) and node != potential_parent:
                hierarchy[node].append(potential_parent)

    # Initialize the PyVis network
    net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "edges": {{
        "smooth": false,
        "color": "gray"
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{
          "enabled": true
        }}
      }},
      "layout": {{
        "hierarchical": {{
          "enabled": false,
          "direction": "DU",
          "sortMethod": "directed"
        }}
      }}
    }}
    """)

    # Add nodes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.coolwarm(1 - normalized_time)) if not pd.isna(normalized_time) else default_color
        net.add_node(
            str(node),
            label=f"{node}",
            title=f"Freq: {freq}, Avg Time: {weighted_time:.2f}" if not pd.isna(weighted_time) else "N/A",
            color=color,
            size=freq * node_size_scale
        )

    # # Add edges based on hierarchy
    # for node, parents in hierarchy.items():
    #     for parent in parents:
    #         net.add_edge(str(parent), str(node), width=1, color="black", title="Subset Relationship")

    # Add transitions as additional edges
    for (source, target), weight in transitions_freq.items():
        net.add_edge(str(source), str(target), weight=weight * edge_width_scale, title=f"Freq: {weight}", color="gray")

    # Save the visualization to an HTML file 
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename

########################################################################################
########################################################################################
########################################################################################
# # # Original: plot_abstraction_graph_movie
def plot_abstraction_graph_movie(abstract_time_series, interval=400, title="Abstraction Graph - Time"):
    """
    Creates an animated abstraction graph that evolves over time.

    The animation highlights:
    - The current node in red.
    - Past nodes gradually fading (transparent) as time progresses.
    - The last three unique abstraction sets in varying sizes and shades of blue.
    - Node size: Depends on the frequency of occurrences of each abstraction.
    - Node color: Represents the weighted average time of occurrence (normalized).
    - Edge thickness: Indicates the frequency of transitions between abstractions.

    Parameters:
    abstract_time_series (pd.DataFrame): A DataFrame with the following columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    interval (int): Interval in milliseconds between frames of the animation. Defaults to 400.
    title (str): Title of the graph animation. Defaults to "Abstraction Graph - Time".

    Returns:
    None. Displays an animation of the abstraction graph.
    """
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last')  # Keep last occurrences of unique abstractions
    last_three_unique = unique_abstractions.tail(3)  # Last three unique sets

    nodes = abstractions.apply(tuple).value_counts(sort=False)
    # Ordenar el Series por la longitud de la primera coordenada de las claves
    # nodes= nodes.sort_index(key=lambda x: [len(k) if isinstance(k, (list, tuple)) else 0 for k in x])

    print(nodes)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.Blues(normalized_time))  # Convert the color to hexadecimal
        G.add_node(node, size=freq, color=color)

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Hierarchical layout
    pos = nx.planar_layout(G) # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    pos = {node: (x, -y) for node, (x, y) in pos.items()}

    # Create the figure and axis for animation
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis("off")

    # Sizes and colors for the last three nodes
    highlight_sizes = [1 / 4, 1 / 2, 1]
    highlight_colors = ['lightblue', 'dodgerblue', 'blue']  # Different shades of blue

    # Create a function for drawing each frame
    def update_frame(i):
        ax.clear()  # Clear the current frame

        # Draw the nodes
        current_abstraction = abstractions.iloc[i]
        visited_nodes = abstractions.iloc[:i + 1].apply(tuple).value_counts().index
        node_sizes = [G.nodes[node]["size"] * 100 for node in G.nodes]

        # Adjust node colors based on visited state, with transparency for past nodes
        node_colors = [
            G.nodes[node]["color"] if node not in visited_nodes else (0.8, 0.8, 1, 0.3)  # Fade visited nodes
            for node in G.nodes
        ]

        # Highlight the last three unique nodes
        for idx, last_node in enumerate(last_three_unique):
            if last_node in G.nodes:
                node_index = list(G.nodes).index(last_node)
                node_sizes[node_index] = highlight_sizes[idx] * 500  # Update size
                node_colors[node_index] = highlight_colors[idx]  # Update color

        # Draw the nodes
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, ax=ax)

        # Draw the edges with width based on transition frequency
        edge_weights = [G.edges[edge]["weight"] for edge in G.edges]
        nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.7, ax=ax)

        # Draw the current node as a larger, fully opaque red node
        current_node = tuple(current_abstraction)
        nx.draw_networkx_nodes(G, pos, nodelist=[current_node], node_size=500, node_color='red', ax=ax)

        # Add labels for nodes
        node_labels = {node: f"({', '.join(map(str, node))})" for node in G.nodes}
        nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_color="black", ax=ax)

        # Title
        ax.set_title(f"{title} {i}", fontsize=15)

        return ax,

    # Create the animation by updating the frame for each timestep
    ani = animation.FuncAnimation(fig, update_frame, frames=len(abstractions), interval=interval, blit=False, repeat=False)

    # Show the animation
    plt.show()

##############################################################################################

def plot_abstraction_graph_movie_3last_nodes(abstract_time_series, 
                                             first_node_size=125, opacity_first_node=0.5,
                                             last_node_size=500, last_adjust_sizes = [1/4, 1/2, 1], 
                                             colour_last_nodes = ['blue', 'blue', 'blue'], 
                                             interval=400, 
                                             title="Abstraction Graph - Time"):
    """
    Creates an animated abstraction graph showing the evolution of abstraction sets over time.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame with the following columns:
        - 'Time': Numeric, indicating the occurrence time of each abstraction.
        - 'Abstraction': Iterable, representing abstraction sets.
    
    first_node_size : int, optional
        Size of the nodes that are not part of the last three unique abstractions. Defaults to 125.
    
    opacity_first_node : float, optional
        Opacity of nodes representing visited states in the animation. Range [0, 1]. Defaults to 0.5.
    
    last_node_size : int, optional
        Base size for the last three unique nodes. Their sizes are scaled by `last_adjust_sizes`. Defaults to 500.
    
    last_adjust_sizes : list of float, optional
        Scaling factors for the sizes of the last three unique nodes. Defaults to [1/4, 1/2, 1].
    
    colour_last_nodes : list of str, optional
        Colors for the last three unique nodes. Should contain three color specifications. Defaults to ['blue', 'blue', 'blue'].
    
    interval : int, optional
        Interval in milliseconds between frames of the animation. Defaults to 400.
    
    title : str, optional
        Title of the graph animation. Defaults to "Abstraction Graph - Time".

    Returns:
    -------
    None
        Displays an animation of the abstraction graph.
    """

    abstractions = abstract_time_series["Abstraction"]
    abstractions= pd.concat([abstractions, pd.Series([abstractions.iloc[-1]])], ignore_index=True)
    times = abstract_time_series["Time"]
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last')  # Keep last occurrences of unique abstractions
    last_three_unique = unique_abstractions.tail(3)  # Last three unique sets
    nodes = abstractions.apply(tuple).value_counts(sort=False) # Count node frequencies and without sort (sort=False).                 # nodes = nodes.sort_index(key=lambda x: [len(k) if isinstance(k, (list, tuple)) else 0 for k in x])
    
    # Create unique labels for nodes (S1, S2, S3, ...)
    node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}
    # print(node_labels)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        nodee = ', '.join(node)                    # Delete inverted commas from the elements of the vector nodes
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.Blues(normalized_time))       # Convert the color to hexadecimal
        label = f"{node_labels[node]} ({freq}) \n[{nodee}]" # Node labels: S1 (2) [a, b,c]
        G.add_node(node, size=freq, color=color, label=label)

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Hierarchical layout
    pos = nx.planar_layout(G) # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    pos = {node: (x, -y) for node, (x, y) in pos.items()}

    # Create the figure and axis for animation
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.axis("off")

    # Create a function for drawing each frame
    def update_frame(i):
        ax.clear()     # Limpia el eje para la nueva iteración
        ax.axis("off") # Desactiva el borde y las marcas del eje

        # Draw the nodes
        current_abstraction = abstractions.iloc[i]
        visited_nodes = abstractions.iloc[:i + 1].apply(tuple).value_counts().index
        node_sizes = [G.nodes[node]["size"] * first_node_size for node in G.nodes]

        # Adjust node colors based on visited state, with transparency for past nodes
        node_colors = [
            G.nodes[node]["color"] if node not in visited_nodes else (0.8, 0.8, 1, opacity_first_node)  # Fade visited nodes. Azul claro:0.8, 0.8, 1. Opacidad: 0.5
            for node in G.nodes
        ]

        for idx, last_node in enumerate(last_three_unique):
            if last_node in G.nodes:
                node_index = list(G.nodes).index(last_node)
                node_sizes[node_index] = last_adjust_sizes[idx] * last_node_size  # Adjust size
                node_colors[node_index] = colour_last_nodes[idx]      # Uniform blue color

        # Draw the nodes
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, ax=ax)

        # Draw the edges with width based on transition frequency
        edge_weights = [G.edges[edge]["weight"] for edge in G.edges]
        nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.9, ax=ax)

        # Draw the current node as a larger, fully opaque red node, except for the last frame
        if i < len(abstractions)-1:  # Evita el último frame
            current_node = tuple(current_abstraction)
            nx.draw_networkx_nodes(G, pos, nodelist=[current_node], node_size=500, node_color='red', ax=ax)
        
        # Add labels for nodes 
        label_pos = {node: (x, y + 10) for node, (x, y) in pos.items()}  # Ajusta '10' según el desplazamiento deseado        
        nx.draw_networkx_labels(
            G, label_pos, labels=nx.get_node_attributes(G, 'label'), 
            font_size=8, font_color="black", ax=ax, 
            bbox=dict(facecolor="white", edgecolor="white", boxstyle="round,pad=0.3", alpha=0)
        )
        
        # Title
        if i < len(abstractions) - 1:  # Para todos los frames excepto el último
            ax.set_title(f"{title} {i}", fontsize=15)
        else:  # Para el último frame
            j=len(abstractions)-2
            ax.set_title(f"{title} {j}", fontsize=15)

        return ax,

    # Create the animation by updating the frame for each timestep
    ani = animation.FuncAnimation(fig, update_frame, frames=len(abstractions), interval=interval, blit=False, repeat=False)

    # Show the animation 
    plt.show()

##############################################################################################    

def get_plot_abstraction_graph_movie_html(abstract_time_series, 
                                          first_node_size=25, colour_first_nodes="lightcyan",
                                          last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], 
                                          colour_last_nodes="red", 
                                          interval=1000, showlegend=True,
                                          title="Abstraction Graph Movie",
                                          filename="abstraction_graph_movie.html"):
    """
    Generates an animated HTML file visualizing the abstraction graph over time.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Time": Corresponding timestamps for each abstraction state.        
        - "Abstraction": A list of abstraction states over time.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_first_nodes : str, optional. Color for the nodes, by default "lightcyan".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : list of str, optional. Colors for the last three unique nodes, by default ['blue', 'blue', 'blue'].
    interval : int, optional. Time interval between frames in milliseconds, by default 1000.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Abstraction Graph Movie".
    filename : str, optional. Name of the output HTML file, by default "abstraction_graph_movie.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
    """
    abstractions = abstract_time_series["Abstraction"]                           # Obtiene la columna "Abstraction" del DataFrame 'abstract_time_series'.
    abstractions = pd.concat([abstractions, pd.Series([abstractions.iloc[-1]])], ignore_index=True)  # Añade el último valor de la columna al final, para garantizar que se considere la transición final.
    times = abstract_time_series["Time"]                                         # Obtiene la columna "Time" del DataFrame 'abstract_time_series'.
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last') # Convierte los elementos de 'abstractions' a tuplas y elimina duplicados, manteniendo el último.
    last_three_unique = unique_abstractions.tail(3)                              # Selecciona las tres últimas abstracciones únicas.
    nodes = abstractions.apply(tuple).value_counts(sort=False)                   # Cuenta las ocurrencias de cada tupla única en 'abstractions'.

    # Create unique labels for nodes
    node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}  # Crea etiquetas únicas para los nodos, como "S1", "S2", etc.

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  # Crea las transiciones entre los elementos consecutivos de 'abstractions'.
    transitions_freq = pd.Series(transitions).value_counts()  # Cuenta la frecuencia de cada transición.

    # Create the graph
    G = nx.DiGraph()  # Crea un grafo dirigido vacío.

    for node, freq in nodes.items():                          # Itera sobre los nodos y sus frecuencias.
        nodee = ', '.join(node)                               # Convierte la tupla del nodo en una cadena de texto.
        color = colour_first_nodes                            # Define el color de los nodos iniciales.
        label = f"{node_labels[node]} ({freq})"               # Etiqueta del nodo con su frecuencia.
        hover_info = f"[{nodee}]"                             # Información detallada para el hover.
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)  # Añade el nodo al grafo con atributos adicionales.

    for (source, target), weight in transitions_freq.items(): # Itera sobre las transiciones y sus frecuencias.
        G.add_edge(source, target, weight=weight)             # Añade las transiciones como aristas con su peso (frecuencia).

    # Primero necesitas asignar niveles a los nodos basado en el tamaño del conjunto
    # En la función get_plot_abstraction_graph_movie_html, alrededor de línea 1554:

    # Crear el subset_key correctamente para multipartite_layout
    subset_key = {}
    max_size = 0

    # Primero encontrar el tamaño máximo entre todos los nodos
    for node in G.nodes:
        if isinstance(node, (tuple, list, set)):
            node_size = len(node)
            max_size = max(max_size, node_size)
        else:
            max_size = max(max_size, 1)

    # Asignar niveles (conjuntos más grandes nivel más alto)
    for node in G.nodes:
        if isinstance(node, (tuple, list, set)):
            node_size = len(node)
            level = node_size  # Conjuntos más grandes = nivel más alto
        else:
            level = 1  # Nodos individuales
        
        # Agrupar nodos por nivel
        if level not in subset_key:
            subset_key[level] = []
        subset_key[level].append(node)

    # Usar el layout multipartita
    pos = nx.multipartite_layout(G, subset_key=subset_key, align='vertical')  # Cambiar a vertical

    # INVERTIR LAS COORDENADAS para que los conjuntos grandes estén arriba
    for node in pos:
        x, y = pos[node]
        # Intercambiar X e Y para orientación vertical
        # e invertir Y para que niveles más altos estén arriba
        pos[node] = (-y, x)  # Intercambiar y ajustar orientación

    # Ajustar espaciado para mejor visualización
    for node in pos:
        x, y = pos[node]
        pos[node] = (x * 1.2, y * 1.5)  # Aumentar separación
    # # pos = nx.planar_layout(G) # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot") # Calcula las posiciones de los nodos utilizando Graphviz.
    # pos = {node: (x, -y) for node, (x, y) in pos.items()}     # Invierte la coordenada y para ajustar la visualización.

    frames = []                                               # Lista para almacenar los frames de la animación.
    sliders_steps = []                                        # Lista para los pasos del control deslizante de la animación.
    
    for i in range(len(abstractions)-1):                      # Itera sobre las abstracciones (menos el último valor).
        node_sizes = [G.nodes[node]["size"] * first_node_size for node in G.nodes]  # Define el tamaño de los nodos.
        node_colors = [G.nodes[node]["color"] for node in G.nodes]                  # Define el color de los nodos.

        # # Destacar los nodos actuales y ajustar tamaños progresivamente
        # for offset, size_factor in enumerate(last_adjust_sizes[::-1]):  # Ajustes de tamaño
        #     # Calcular el índice máximo de los nodos a graficar
        #     max_nodes_to_plot = i + 1  # Número de nodos a graficar en función del tiempo `i`
            
        #     # Verificar que el nodo actual esté dentro del rango permitido
        #     if offset < max_nodes_to_plot:
        #         current_node = tuple(abstractions.iloc[i - offset])  # Obtener el nodo actual
                
        #         # Verificar si el nodo existe en el grafo
        #         if current_node in G.nodes:
        #             index = list(G.nodes).index(current_node)        # Índice del nodo
        #             node_sizes[index] = last_node_size * size_factor # Ajustar tamaño
        #             node_colors[index] = colour_last_nodes           # Cambiar color del nodo
        #         else:
        #             continue  # Saltar si el nodo no existe en el grafo
        # # Destacar los nodos actuales y ajustar tamaños progresivamente
        # for offset, size_factor in enumerate(last_adjust_sizes[::-1]):  # Ajustes de tamaño
        #     # Calcular el índice máximo de los nodos a graficar
        #     max_nodes_to_plot = i + 1  # Número de nodos a graficar en función del tiempo `i`
            
        #     # Verificar que el nodo actual esté dentro del rango permitido
        #     if offset < max_nodes_to_plot:
        #         current_node = tuple(abstractions.iloc[i - offset])  # Obtener el nodo actual
                
        #         # Verificar si el nodo existe en el grafo
        #         if current_node in G.nodes:
        #             index = list(G.nodes).index(current_node)  # Índice del nodo
                    
        #             # Encontrar el tamaño máximo de nodos repetidos
        #             max_size = max(node_sizes[index], last_node_size * size_factor)
        #             node_sizes[index] = max_size  # Ajustar tamaño al máximo
        #             node_colors[index] = colour_last_nodes  # Cambiar color del nodo
        #         else:
        #             continue  # Saltar si el nodo no existe en el grafo
        # Destacar los nodos actuales y ajustar tamaños progresivamente
        for offset, size_factor in enumerate(last_adjust_sizes[::-1]):  # Ajustes de tamaño
            # Calcular el índice máximo de los nodos a graficar
            max_nodes_to_plot = i + 1  # Número de nodos a graficar en función del tiempo `i`
            
            # Verificar que el nodo actual esté dentro del rango permitido
            if offset < max_nodes_to_plot:
                current_node = tuple(abstractions.iloc[i - offset])  # Obtener el nodo actual
                
                # Verificar si el nodo existe en el grafo
                if current_node in G.nodes:
                    index = list(G.nodes).index(current_node)  # Índice del nodo
                    
                    # Condición para tres nodos consecutivos diferentes
                    if offset >= 2:  # Verificar si hay al menos tres nodos consecutivos
                        prev_node_1 = tuple(abstractions.iloc[i - offset + 1])
                        prev_node_2 = tuple(abstractions.iloc[i - offset + 2])
                        if (
                            prev_node_1 in G.nodes
                            and prev_node_2 in G.nodes
                            and current_node != prev_node_1
                            and current_node != prev_node_2
                            and prev_node_1 != prev_node_2
                        ):
                            node_sizes[index] = last_node_size * size_factor
                        else:
                            # Ajustar tamaño al máximo existente
                            max_size = max(node_sizes[index], last_node_size * size_factor)
                            node_sizes[index] = max_size
                    else:
                        # Ajustar tamaño al máximo existente
                        max_size = max(node_sizes[index], last_node_size * size_factor)
                        node_sizes[index] = max_size
                    
                    node_colors[index] = colour_last_nodes  # Cambiar color del nodo
                else:
                    continue  # Saltar si el nodo no existe en el grafo

 
        edge_x = []                    # Lista para las coordenadas x de las aristas.
        edge_y = []                    # Lista para las coordenadas y de las aristas.
        annotations = []               # Lista para las anotaciones de las aristas.

        for source, target in G.edges: # Itera sobre las aristas.
            x0, y0 = pos[source]       # Obtiene las coordenadas del nodo fuente.
            x1, y1 = pos[target]       # Obtiene las coordenadas del nodo objetivo.
            edge_x += [x0, x1, None]   # Añade las coordenadas de la arista.
            edge_y += [y0, y1, None]   # Añade las coordenadas de la arista.

            annotations.append(dict(   # Añade anotaciones para cada arista.
                ax=x0, ay=y0,
                x=x1, y=y1,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=2,
                arrowwidth=1,
                arrowcolor="gray"
            ))

        frames.append(go.Frame(                                  # Añade un nuevo frame para la animación.
            data=[
                go.Scatter(x=[pos[node][0] for node in G.nodes], # Posiciones x de los nodos.
                        y=[pos[node][1] for node in G.nodes],    # Posiciones y de los nodos.
                        mode='markers+text',                     # Modo de visualización de los nodos.
                        marker=dict(size=node_sizes, color=node_colors, opacity=1), # Tamaño y color de los nodos.
                        text=[G.nodes[node]["label"] for node in G.nodes],          # Etiquetas de los nodos.
                        textposition="bottom center"),           # Posición del texto.
                go.Scatter(x=edge_x, y=edge_y,                   # Coordenadas de las aristas.
                        mode='lines',                            # Modo de visualización de las aristas.
                        line=dict(width=1, color='gray'))        # Estilo de las aristas.
            ],
            layout=dict(annotations=annotations),                # Añade las anotaciones a la visualización.
            name=f"Frame {i}"                                    # Nombre del frame para la animación.
        ))
        sliders_steps.append(dict(                               # Añade un paso al control deslizante de la animación.
            args=[[f"Frame {i}"], dict(frame=dict(duration=0, redraw=True), mode="immediate")],  # Define el paso.
            label=f"{i}",                                        # Etiqueta del paso.
            method="animate"                                     # Define el método de animación.
        ))

    legend_text = "<br>".join([f"{node_labels[node]}= [{', '.join(node)}]" for node in G.nodes])  # Genera el texto de la leyenda con etiquetas de nodos (S1, S2, ...) y sus elementos.

    fig = go.Figure(
        data=[
            go.Scatter(
                x=[pos[node][0] for node in G.nodes],  # Coordenadas X de los nodos.
                y=[pos[node][1] for node in G.nodes],  # Coordenadas Y de los nodos.
                mode='markers+text',                   # Modo de visualización: muestra puntos (markers) y texto (text).
                marker=dict(size=node_sizes, color=node_colors, opacity=1),  # Configuración de los nodos: tamaño, color, y opacidad.
                text=[G.nodes[node]["label"] for node in G.nodes],  # Etiquetas que se muestran en cada nodo.
                textposition="bottom center",          # Posición del texto respecto al nodo (debajo y centrado).
                hovertext=[G.nodes[node]["hovertext"] for node in G.nodes],  # Texto que aparece al pasar el mouse sobre un nodo.
                hoverinfo="text",                      # Especifica que se muestra el texto de hover al interactuar.
                name=legend_text                       # Asocia el texto de la leyenda con el gráfico.
            )
        ],
        layout=go.Layout(
            title=title,                      # Título del gráfico.
            updatemenus=[                     # Configuración de los botones de control (play y pause).
                dict(
                    type="buttons",           # Define un grupo de botones.
                    showactive=False,         # Desactiva el resaltado del botón seleccionado.
                    buttons=[
                        dict(
                            label="Play",     # Etiqueta del botón "Play".
                            method="animate", # Método que activa la animación.
                            args=[None, dict(frame=dict(duration=interval, redraw=True), fromcurrent=True)]  # Configuración de la animación.
                        ),
                        dict(
                            label="Pause",    # Etiqueta del botón "Pause".
                            method="animate", # Método que pausa la animación.
                            args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")]  # Configuración de pausa.
                        )
                    ]
                )
            ],
            sliders=[{  # Configuración del control deslizante para avanzar manualmente entre los cuadros de la animación.
                'steps': sliders_steps,   # Pasos del deslizador.
                'currentvalue': {
                    'prefix': "Time: ",   # Prefijo que aparece junto al valor actual del deslizador.
                    'font': {'size': 16}, # Tamaño de fuente del prefijo.
                    'visible': True,      # Muestra el valor actual.
                },
                'x': 0.1,   # Posición horizontal del deslizador.
                'len': 0.9, # Longitud del deslizador en proporción al gráfico.
            }],
            xaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje X.
            yaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje Y.
            plot_bgcolor='rgba(0,0,0,0)', # Fondo transparente del gráfico.
            width=1400,                   # Ancho del gráfico.
            height=800,                   # Altura del gráfico.
            annotations=annotations,      # Anotaciones adicionales para destacar elementos.
            showlegend=showlegend,        # Habilita la visualización de la leyenda.
            legend=dict( 
                title=dict(text="Node Legend", font=dict(family="Arial", size=14, color="darkblue")),  # Título de la leyenda.
                bgcolor="white",  # Fondo blanco.
                bordercolor="white",                               # Borde blanco.
                borderwidth=2,                                     # Grosor del borde.
                font=dict(family="Arial", size=12, color="black"), # Fuente de los elementos.
                x=1,                                               # Posición horizontal a la derecha.
                y=1.0                                              # Posición vertical en la parte superior.
            )
        ),
        frames=frames  # Cuadros de la animación.
    )

    fig.write_html(filename)
    return filename


############################################################################################## 
 
import os
import webbrowser

def plot_abstraction_graph_movie_html(
                abstract_time_series, 
                first_node_size=25, colour_first_nodes="lightcyan",
                last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], 
                colour_last_nodes="red", 
                interval=1000, showlegend=True,
                title="Abstraction Graph Movie",
                filename="abstraction_graph_movie.html"):
    """
    Generates an interactive animation of an abstraction graph from an abstract time series, 
    saving the animation as an HTML file in the specified directory. The animation shows how 
    the graph evolves over time, highlighting the most important nodes.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Abstraction": A list of abstraction states over time.
        - "Time": Corresponding timestamps for each abstraction state.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_first_nodes : str, optional. Color for the nodes, by default "lightcyan".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : list of str, optional. Colors for the last three unique nodes, by default ['blue', 'blue', 'blue'].
    interval : int, optional. Time interval between frames in milliseconds, by default 1000.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Abstraction Graph Movie".
    filename : str, optional. Name of the output HTML file, by default "abstraction_graph_movie.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
    Example:
    -------
    # Load an abstraction time series
    abstract_time_series = pd.DataFrame({
        'Time': [0, 1, 2],
        'Abstraction': [['A', 'B'], ['B'], ['C', 'D']]        
    })

    # Call the function to generate the animation
    plot_abstraction_graph_movie_html(abstract_time_series, filename="graph_movie.html")
    """    
    # Create the directory if it doesn't exist
    output_dir = "visualizations/plot_abstraction_graph_movie_html"
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct the full file path
    filepath = os.path.join(output_dir, filename)
    
    # Call the `rn_get_visualization` function to generate the HTML file with the correct path
    get_plot_abstraction_graph_movie_html(
        abstract_time_series,
        first_node_size=first_node_size, colour_first_nodes=colour_first_nodes,
        last_node_size=last_node_size, last_adjust_sizes=last_adjust_sizes, 
        colour_last_nodes=colour_last_nodes, 
        interval=interval, showlegend=showlegend,
        title=title,
        filename=filepath  # Pass the full path here
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filepath) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `rn_get_visualization` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The animation was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")
    
    return abs_path  # Return the absolute path for reference

##############################################################################################    













##############################################################################################

def plot_abstraction_graph_movie_html0(abstract_time_series, interval=400, title="Abstraction Graph - Time", 
                                      filename="abstraction_graph_movie0.html"):
    """
    Creates an interactive abstraction graph that evolves over time and saves it as an HTML file.

    Parameters:
    abstract_time_series (pd.DataFrame): A DataFrame with the following columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    interval (int): Interval in milliseconds between frames of the animation. Defaults to 400.
    title (str): Title of the graph animation. Defaults to "Abstraction Graph - Time".
    filename (str): The name of the HTML file to save the animation.

    Returns:
    None. Saves the interactive abstraction graph animation to an HTML file.
    """
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last')  # Keep last occurrences of unique abstractions
    last_three_unique = unique_abstractions.tail(3)  # Last three unique sets

    nodes = abstractions.apply(tuple).value_counts()

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.Blues(normalized_time))  # Convert the color to hexadecimal
        G.add_node(node, size=freq, color=color)

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Hierarchical layout
    pos = nx.planar_layout(G) # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    pos = {node: (x, -y) for node, (x, y) in pos.items()}

    # Prepare the PyVis network
    net = Network(height="750px", width="100%", directed=True)

    # Add nodes to the PyVis network
    for node in G.nodes:
        node_size = G.nodes[node]["size"] * 100
        color = G.nodes[node]["color"]
        
        # Add node to the network (convert to string for PyVis compatibility)
        net.add_node(str(node), label=f"({', '.join(map(str, node))})", title=f"({', '.join(map(str, node))})", color=color, size=node_size, x=pos[node][0], y=pos[node][1])

    # Add edges to the PyVis network
    for source, target in G.edges:
        weight = G.edges[source, target].get("weight", 1)
        # Add edges to network (convert nodes to string)
        net.add_edge(str(source), str(target), width=weight)

    # Set up the layout options for the network
    net.set_options("""
    {
        "nodes": {
            "shape": "dot",
            "font": {
                "size": 14,
                "align": "center"
            },
            "borderWidth": 2,
            "borderColor": "black"
        },
        "edges": {
            "smooth": false,
            "color": "#888888",
            "width": 1
        },
        "physics": {
            "enabled": false
        },
        "layout": {
            "hierarchical": {
                "enabled": true,
                "direction": "DU",
                "sortMethod": "directed"
            }
        }
    }
    """) 
        # Save the visualization to an HTML file 
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename

##############################################################################################
# # Películas animadas de gráficos de abstracción y semiorganizaciones
##############################################################################################

# # Function to generate a Hasse diagram of the set of species in a reaction network.
def plot_hasse_diagram(species_set, node_size=25, color_node='cyan', arrowcolor='black', bgcolor_legend="white", bordercolor_legend="white", name_legend="Legend:", title="Hasse diagram"):
    """
    Generates and plots a Hasse diagram for a given set of species. The diagram includes node labels 
    and a legend, showing the relationships between different subsets of species. 

    Parameters:
    ----------
    species_set : set
        A set of species that will be used to generate the Hasse diagram.
    node_size : int, optional
        The size of the nodes in the diagram, by default 25.
    color_node : str, optional
        The color for the nodes in the diagram, by default 'cyan'.
    arrowcolor : str, optional
        The color of the arrows connecting the nodes, by default 'black'.
    bgcolor_legend : str, optional
        The background color of the legend, by default 'white'.
    bordercolor_legend : str, optional
        The border color of the legend, by default 'white'.
    name_legend : str, optional
        The title for the legend, by default 'Legend:'.
    title : str, optional
        The title of the diagram, by default 'Hasse diagram'.

    Returns:
    -------
    None
        The function displays the Hasse diagram plot using Plotly.
    
    Example:
    -------
    species_set = {"A", "B", "C", "D"}
    plot_hasse_diagram(species_set)
    """
    # Generar todos los subconjuntos del conjunto dado
    elementos = list(species_set)
    subconjuntos = [set(comb) for i in range(len(elementos) + 1) for comb in combinations(elementos, i)]
    
    # Crear nodos y relaciones
    nodos = []
    edges = []
    niveles = {}

    for subconjunto in subconjuntos:
        nodo = frozenset(subconjunto)
        nodos.append(nodo)
        niveles[nodo] = len(subconjunto)
    
    for subconjunto1 in subconjuntos:
        for subconjunto2 in subconjuntos:
            if subconjunto1 < subconjunto2 and len(subconjunto2 - subconjunto1) == 1:
                edges.append((frozenset(subconjunto1), frozenset(subconjunto2)))

    # Etiquetar nodos
    etiquetas = {nodo: f"S{idx + 1}" for idx, nodo in enumerate(nodos)}
    
    # Posiciones de los nodos (usando el nivel como coordenada y distribuyendo en X)
    posiciones = {}
    nivel_max = max(niveles.values())
    for nivel in range(nivel_max + 1):
        nodos_nivel = [nodo for nodo in nodos if niveles[nodo] == nivel]
        x_coords = range(-len(nodos_nivel) // 2, len(nodos_nivel) // 2 + 1)
        for x, nodo in zip(x_coords, nodos_nivel):
            posiciones[nodo] = (x, nivel)

    # Datos para Plotly
    x_coords = []
    y_coords = []
    labels = []
    hovertexts = []
    for nodo, (x, y) in posiciones.items():
        x_coords.append(x)
        y_coords.append(y)  # Sin invertir Y, el vacío queda abajo
        labels.append(etiquetas[nodo])
        hovertexts.append(f"{[', '.join(nodo) if nodo else '∅']}")

    # Crear las aristas con flechas
    annotations = []
    for nodo1, nodo2 in edges:
        x0, y0 = posiciones[nodo1]
        x1, y1 = posiciones[nodo2]
        annotations.append(dict(
            ax=x0, ay=y0,
            x=x1, y=y1,
            xref='x', yref='y',
            axref='x', ayref='y',
            showarrow=True,
            arrowhead=2,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor=arrowcolor
        ))

    # Trazar los nodos
    node_trace = go.Scatter(
        x=x_coords, y=y_coords,
        mode='markers+text',
        text=labels,
        hovertext=hovertexts,
        textposition="top center",
        marker=dict(
            size=node_size,
            color=color_node,
            line=dict(width=2, color=color_node)
        )
    )

    # Crear la leyenda
    leyenda_texto = "<br>".join([f"{etiquetas[nodo]}= [{', '.join(nodo) if nodo else '∅'}]" for nodo in nodos])

    # Crear la figura
    fig = go.Figure(data=[node_trace])
    fig.update_layout(
        title=title,
        showlegend=False,
        xaxis=dict(
            showgrid=False, zeroline=False, showticklabels=False
        ),
        yaxis=dict(
            showgrid=False, zeroline=False, showticklabels=False
        ),
        plot_bgcolor='white',
        annotations=annotations
    )
    
    # Agregar la leyenda como anotación
    fig.add_annotation(
        x=1.05, y=1, xref="paper", yref="paper",
        text=f"<b>{name_legend}</b><br>{leyenda_texto}",
        showarrow=False,
        align="left",
        bordercolor=bordercolor_legend,
        borderwidth=1,
        bgcolor=bgcolor_legend, 
        font=dict(size=12)
    )

    fig.show()

##############################################################################################
# # Function to generate an interactive graph from a Hasse diagram of Semi-Organisations and Abstractions.
def get_film_semiorganizations_abstractions_html(abstract_time_series, input_data,
                                          first_node_size=25, 
                                          colour_abs_nodes="lightcyan", colour_semiorg_nodes="lightgreen", colour_cap_nodes="orange",
                                          last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], colour_last_nodes="red", 
                                          interval=100, showlegend=True,
                                          title="Semi-organisations and Abstractions Film",
                                          filename="semiorg_abstraction_film.html"):
    """
    Generates an animated HTML file visualizing the abstraction graph over time.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Time": Corresponding timestamps for each abstraction state.        
        - "Abstraction": A list of abstraction states over time.
    input_data : list of lists
        A list of subsets to visualize as additional nodes and intersections.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_abs_nodes : str, optional. Color for the nodes, by default "lightcyan".
    colour_semiorg_nodes : str, optional. Color for semi-organization nodes, by default "lightgreen".
    colour_cap_nodes : str, optional. Color for the intersection nodes, by default "orange".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : str, optional. Colors for the last three unique nodes, by default ['red', 'red', 'red'].
    interval : int, optional. Time interval between frames in milliseconds, by default 100.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Semi-organisations and Abstractions Film".
    filename : str, optional. Name of the output HTML file, by default "semiorg_abstraction_film.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
    """
    # Convert input_data to unique sets
    unique_subsets = []
    for sublist in input_data:
        if set(sublist) not in [set(x) for x in unique_subsets]:
            unique_subsets.append(sublist) 

    # Create a list of unique sets for reference
    set_of_sets = [tuple(subset) for subset in unique_subsets] 
    # Ordenar set_of_sets por la cantidad de elementos en cada subconjunto
    set_of_sets = sorted([tuple(sorted(subset)) for subset in set_of_sets], key=len) 

    set_of_sets2 = pd.Series(set_of_sets) # Convertir a Series para aplicar .apply() y .value_counts()
    set_of_sets2 = set_of_sets2.apply(tuple).value_counts(sort=False) # Luego puedes aplicar .apply() y .value_counts()
    # print("set_of_sets2")
    # print(set_of_sets2)
    # Create set names based on their level
    Set_names = {node: f"S{idx + 1}" for idx, node in enumerate(set_of_sets2.keys())}
    # Set_names = [f"S{i+1}" for i in range(len(set_of_sets2))]

    # Create set labels to display when hovering over nodes
    abstractions_0 = abstract_time_series["Abstraction"] # Abstraction states over time # abstractions = pd.concat([abstractions, pd.Series([abstractions.iloc[-1]])], ignore_index=True) # Add the last state to the end
    abstractions_0 = abstractions_0.apply(sorted)
    # print("abstractions_0")
    # print(abstractions_0)
    
    abstractions = pd.concat([abstractions_0, pd.Series(set_of_sets)], ignore_index=True)
    # unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last') # Unique abstractions times = abstract_time_series["Time"] # Time values for each abstraction state # last_three_unique = unique_abstractions.tail(3) # Get the last three unique abstractions

    # Ordenar los elementos de cada subconjunto en unique_subsets
    unique_subsets = [sorted(subset) for subset in unique_subsets]
    # print()
    # print("unique_abstractions", unique_abstractions)

    nodes = abstractions.apply(tuple).value_counts(sort=False) 
    nodes = nodes[~nodes.index.duplicated()]

    nodes_0 = abstractions_0.apply(tuple).value_counts(sort=False)
    nodes_1 = nodes.loc[~nodes.index.isin(nodes_0.index)] 

    # print("nodes")
    # print(nodes)  
    # print("nodes_0")  
    # print(nodes_0)
    # print("nodes_1")    
    # print(nodes_1)     

    # Create unique labels for nodes
    # node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}
    
    # Crear node_labels para nodes_0 y nodes_1
    node_labels_0 = {node: f"A{idx + 1}" for idx, node in enumerate(nodes_0.keys())}
    node_labels_1 = {node: f"S{idx + 1}" for idx, node in enumerate(nodes_1.keys())}

    # Fusionar los diccionarios
    node_labels = {**node_labels_0, **node_labels_1}
    # print()
    # print(node_labels)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    for node, freq in nodes_0.items():
        nodee = ', '.join(node)
        color = colour_abs_nodes
        label = f"{node_labels_0[node]} ({freq})"
        hover_info = f"[{nodee}]"
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)

    for node, freq in set_of_sets2.items():
        nodee = ', '.join(node)
        color = colour_semiorg_nodes
        label = f"{Set_names[node]} ({freq})"
        hover_info = f"[{nodee}]"
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)

    # Graficar nodos de intersección entre nodes_0 y set_of_sets2
    intersection_nodes = set(nodes_0.keys()).intersection(set_of_sets2.keys())
    for node in intersection_nodes:
        freq_nodes_0 = nodes_0[node]
        # freq_set_of_sets2 = set_of_sets2[node]
        nodee = ', '.join(node)
        color =  colour_cap_nodes
        label = f"{node_labels_0[node]} ∩ {Set_names[node]} ({freq_nodes_0})" #, {freq_set_of_sets2})"
        hover_info = f"[{nodee}]"
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)

    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    # pos = {node: (x, -y) for node, (x, y) in pos.items()}


    # Usar disposición con Graphviz y ajustar con un desplazamiento# pos = nx.shell_layout(G)
    pos = nx.planar_layout(G) # pos = graphviz_layout(G, prog="dot")
    # # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    offset_x, offset_y = 10, 10
    pos = {node: (x + offset_x, y + offset_y) for node, (x, y) in pos.items()}


    frames = []
    sliders_steps = []

    for i in range(len(abstract_time_series["Abstraction"])):
        node_sizes = [G.nodes[node]["size"] * first_node_size for node in G.nodes]
        node_colors = [G.nodes[node]["color"] for node in G.nodes]

        # Destacar los nodos actuales y ajustar tamaños progresivamente
        for offset, size_factor in enumerate(last_adjust_sizes[::-1]):
            max_nodes_to_plot = i + 1

            if offset < max_nodes_to_plot:
                current_node = tuple(abstractions.iloc[i - offset])

                if current_node in G.nodes:
                    index = list(G.nodes).index(current_node)

                    if offset >= 2:
                        prev_node_1 = tuple(abstractions.iloc[i - offset + 1])
                        prev_node_2 = tuple(abstractions.iloc[i - offset + 2])
                        if (
                            prev_node_1 in G.nodes
                            and prev_node_2 in G.nodes
                            and current_node != prev_node_1
                            and current_node != prev_node_2
                            and prev_node_1 != prev_node_2
                        ):
                            node_sizes[index] = last_node_size * size_factor
                        else:
                            max_size = max(node_sizes[index], last_node_size * size_factor)
                            node_sizes[index] = max_size
                    else:
                        max_size = max(node_sizes[index], last_node_size * size_factor)
                        node_sizes[index] = max_size

                    node_colors[index] = colour_last_nodes
                else:
                    continue
 
        edge_x = []                    # Lista para las coordenadas x de las aristas.
        edge_y = []                    # Lista para las coordenadas y de las aristas.
        annotations = []               # Lista para las anotaciones de las aristas.

        for source, target in G.edges: # Itera sobre las aristas.
            x0, y0 = pos[source]       # Obtiene las coordenadas del nodo fuente.
            x1, y1 = pos[target]       # Obtiene las coordenadas del nodo objetivo.
            edge_x += [x0, x1, None]   # Añade las coordenadas de la arista.
            edge_y += [y0, y1, None]   # Añade las coordenadas de la arista.

            annotations.append(dict(   # Añade anotaciones para cada arista.
                ax=x0, ay=y0,
                x=x1, y=y1,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=2,
                arrowwidth=1,
                arrowcolor="gray"
            ))

        frames.append(go.Frame(                                  # Añade un nuevo frame para la animación.
            data=[
                go.Scatter(x=[pos[node][0] for node in G.nodes], # Posiciones x de los nodos.
                        y=[pos[node][1] for node in G.nodes],    # Posiciones y de los nodos.
                        mode='markers+text',                     # Modo de visualización de los nodos.
                        marker=dict(size=node_sizes, color=node_colors, opacity=1), # Tamaño y color de los nodos.
                        text=[G.nodes[node]["label"] for node in G.nodes],          # Etiquetas de los nodos.
                        textposition="bottom center"),           # Posición del texto.
                go.Scatter(x=edge_x, y=edge_y,                   # Coordenadas de las aristas.
                        mode='lines',                            # Modo de visualización de las aristas.
                        line=dict(width=1, color='gray'))        # Estilo de las aristas.
            ],
            layout=dict(annotations=annotations),                # Añade las anotaciones a la visualización.
            name=f"Frame {i}"                                    # Nombre del frame para la animación.
        ))
        sliders_steps.append(dict(                               # Añade un paso al control deslizante de la animación.
            args=[[f"Frame {i}"], dict(frame=dict(duration=0, redraw=True), mode="immediate")],  # Define el paso.
            label=f"{i}",                                        # Etiqueta del paso.
            method="animate"                                     # Define el método de animación.
        ))

    legend_text = "<br>".join([f"{node_labels[node]}= [{', '.join(node)}]" for node in G.nodes])  # Genera el texto de la leyenda con etiquetas de nodos (S1, S2, ...) y sus elementos.
    # legend_text = "<br>".join([f"[{node}]" for node in G.nodes])  # Genera el texto de la leyenda con etiquetas de nodos (S1, S2, ...) y sus elementos.

    fig = go.Figure(
        data=[
            go.Scatter(
                x=[pos[node][0] for node in G.nodes],  # Coordenadas X de los nodos.
                y=[pos[node][1] for node in G.nodes],  # Coordenadas Y de los nodos.
                mode='markers+text',                   # Modo de visualización: muestra puntos (markers) y texto (text).
                marker=dict(size=node_sizes, color=node_colors, opacity=1),  # Configuración de los nodos: tamaño, color, y opacidad.
                text=[G.nodes[node]["label"] for node in G.nodes],  # Etiquetas que se muestran en cada nodo.
                textposition="bottom center",          # Posición del texto respecto al nodo (debajo y centrado).
                hovertext=[G.nodes[node]["hovertext"] for node in G.nodes],  # Texto que aparece al pasar el mouse sobre un nodo.
                hoverinfo="text",                      # Especifica que se muestra el texto de hover al interactuar.
                name=legend_text                       # Asocia el texto de la leyenda con el gráfico.
            )
        ],
        layout=go.Layout(
            title=title,                      # Título del gráfico.
            updatemenus=[                     # Configuración de los botones de control (play y pause).
                dict(
                    type="buttons",           # Define un grupo de botones.
                    showactive=False,         # Desactiva el resaltado del botón seleccionado.
                    buttons=[
                        dict(
                            label="Play",     # Etiqueta del botón "Play".
                            method="animate", # Método que activa la animación.
                            args=[None, dict(frame=dict(duration=interval, redraw=True), fromcurrent=True)]  # Configuración de la animación.
                        ),
                        dict(
                            label="Pause",    # Etiqueta del botón "Pause".
                            method="animate", # Método que pausa la animación.
                            args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")]  # Configuración de pausa.
                        )
                    ]
                )
            ],
            sliders=[{  # Configuración del control deslizante para avanzar manualmente entre los cuadros de la animación.
                'steps': sliders_steps,   # Pasos del deslizador.
                'currentvalue': {
                    'prefix': "Time: ",   # Prefijo que aparece junto al valor actual del deslizador.
                    'font': {'size': 16}, # Tamaño de fuente del prefijo.
                    'visible': True,      # Muestra el valor actual.
                },
                'x': 0.1,   # Posición horizontal del deslizador.
                'len': 0.9, # Longitud del deslizador en proporción al gráfico.
            }],
            xaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje X.
            yaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje Y.
            plot_bgcolor='rgba(0,0,0,0)', # Fondo transparente del gráfico.
            width=1400,                   # Ancho del gráfico.
            height=800,                   # Altura del gráfico.
            annotations=annotations,      # Anotaciones adicionales para destacar elementos.
            showlegend=showlegend,        # Habilita la visualización de la leyenda.
            legend=dict( 
                title=dict(text="Node Legend", font=dict(family="Arial", size=14, color="darkblue")),  # Título de la leyenda.
                bgcolor="white",  # Fondo blanco.
                bordercolor="white",                               # Borde blanco.
                borderwidth=2,                                     # Grosor del borde.
                font=dict(family="Arial", size=12, color="black"), # Fuente de los elementos.
                x=1,                                               # Posición horizontal a la derecha.
                y=1.0                                              # Posición vertical en la parte superior.
            )
        ),
        frames=frames  # Cuadros de la animación.
    )

    fig.write_html(filename)
    return filename

############################################################################################## 
 
def film_semiorganizations_abstractions_html(
                abstract_time_series, input_data,
                first_node_size=25, 
                colour_abs_nodes="lightcyan",  colour_semiorg_nodes="lightgreen", colour_cap_nodes="orange",
                last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], 
                colour_last_nodes="red", 
                interval=100, showlegend=True,
                title="Semi-organisations and Abstractions Film",
                filename="semiorg_abstraction_film.html"):
    """
    Generates an interactive animation of an abstraction graph from an abstract time series, 
    saving the animation as an HTML file. The animation shows how the graph evolves over time, 
    highlighting the most important nodes.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Abstraction": A list of abstraction states over time.
        - "Time": Corresponding timestamps for each abstraction state.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_abs_nodes : str, optional. Color for the nodes, by default "lightcyan".
    colour_semiorg_nodes : str, optional. Color for semi-organization nodes, by default "lightgreen".
    colour_cap_nodes : str, optional. Color for the intersection nodes, by default "orange".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : str, optional. Colors for the last three unique nodes, by default 'red'.
    interval : int, optional. Time interval between frames in milliseconds, by default 100.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Semi-organisations and Abstractions Film".
    filename : str, optional. Name of the output HTML file, by default "semiorg_abstraction_film.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
        
    Example:
    -------
    # Load an abstraction time series
    abstract_time_series = pd.DataFrame({
        'Time': [0, 1, 2],
        'Abstraction': [['A', 'B'], ['B'], ['C', 'D']]        
    })

    # Call the function to generate the animation
    plot_abstraction_graph_movie_html(abstract_time_series, filename="semiorg_abstraction_film.html")
    """
    # Call the `get_film_semiorganizations_abstractions_html` function to generate the HTML file
    get_film_semiorganizations_abstractions_html(
        abstract_time_series, input_data,
        first_node_size=first_node_size, 
        colour_abs_nodes=colour_abs_nodes, colour_semiorg_nodes=colour_semiorg_nodes, colour_cap_nodes=colour_cap_nodes,
        last_node_size=last_node_size, last_adjust_sizes=last_adjust_sizes, 
        colour_last_nodes=colour_last_nodes, 
        interval=interval, showlegend=showlegend,
        title=title,
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    # print("absolute path")
    # print(abs_path)
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `rn_get_visualization` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The animation was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

######################################################################################
# Histograms
######################################################################################

def plot_join_concentration_histogram(time_series, bins=30, alpha=0.7, color='skyblue', edgecolor='black', xlabel="Concentration", ylabel="Frequency", title="Histogram of Species Concentrations"):
    """
    Plots a histogram of combined species concentrations from the time series.

    Parameters:
    ----------
    time_series : pd.DataFrame
        Time series with a 'Time' column and species concentrations as additional columns.
    bins : int, optional
        Number of bins for the histogram (default is 30).
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    color : str, optional
        Color of the histogram bars (default is 'skyblue').
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is 'black').
    xlabel : str, optional
        Label for the x-axis (default is "Concentration").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title : str, optional
        Title of the plot (default is "Histogram of Species Concentrations").
    """
    if 'Time' not in time_series.columns:
        raise ValueError("The DataFrame must include a 'Time' column for time values.")
    
    concentrations = time_series.drop(columns='Time').values.flatten()
    
    plt.figure(figsize=(10, 6))
    plt.hist(concentrations, bins=bins, alpha=alpha, color=color, edgecolor=edgecolor)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_reaction_rate_mak_histogram(reaction_rate_maks, bins=10, color='skyblue', edgecolor='black', alpha=0.7, xlabel="Reaction Rate", ylabel="Frequency", title="Histogram of Reaction Rates"):
    """
    Generates a histogram for the reaction rates of species.

    Parameters:
    ----------
    reaction_rate_maks : pd.DataFrame
        DataFrame with reaction rates for species. Each column represents a species, and rows represent different time points.
    bins : int, optional
        Number of bins in the histogram (default is 10).
    color : str, optional
        Color of the histogram bars (default is 'skyblue').
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is 'black').
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    xlabel : str, optional
        Label for the x-axis (default is "Reaction Rate").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title : str, optional
        Title of the histogram (default is "Histogram of Reaction Rates").
    """
    all_rates = reaction_rate_maks.values.flatten()
    
    plt.figure(figsize=(10, 6))
    plt.hist(all_rates, bins=bins, color=color, edgecolor=edgecolor, alpha=alpha)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_species_histograms(reaction_rate_maks, species_names, bins=10, alpha=0.7, color="skyblue", edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Histogram of"):
    """
    Generates individual histograms for each species in the reaction rates DataFrame.

    Parameters:
    ----------
    reaction_rate_maks : pd.DataFrame
        DataFrame with reaction rates for species.
    species_names : list
        List of species names (column names in the DataFrame).
    bins : int, optional
        Number of bins for the histograms (default is 10).
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    color : str, optional
        Color of the histogram bars (default is "skyblue").
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is "black").
    xlabel : str, optional
        Label for the x-axis (default is "Concentration").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title_plot : str, optional
        Prefix for the titles of each subplot (default is "Histogram of").
    """
    num_species = len(species_names)
    fig, axes = plt.subplots(1, num_species, figsize=(5 * num_species, 4), sharey=True)

    if num_species == 1:
        axes = [axes]

    for ax, species in zip(axes, species_names):
        ax.hist(reaction_rate_maks[species], bins=bins, alpha=alpha, color=color, edgecolor=edgecolor, label=species)
        ax.set_title(f"{title_plot} {species}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()  # Adds the legend to the subplot.

    plt.tight_layout()
    plt.show()


def plot_combined_species_histogram(reaction_rate_maks, species_names, bins=10, alpha=0.7, edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Combined Histogram"):
    """
    Generates a combined histogram of reaction rates for a set of species.

    Parameters:
    ----------
    reaction_rate_maks : pd.DataFrame
        DataFrame with reaction rates for species.
    species_names : list
        List of species names (columns in the DataFrame).
    bins : int, optional
        Number of bins for the histograms (default is 10).
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is "black").
    xlabel : str, optional
        Label for the x-axis (default is "Concentration").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title_plot : str, optional
        Title of the plot (default is "Combined Histogram").
    """
    plt.figure(figsize=(8, 6))

    for species in species_names:
        plt.hist(reaction_rate_maks[species], bins=bins, alpha=alpha, label=species, edgecolor=edgecolor)
    
    plt.title(title_plot)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(title="Species")
    plt.grid(axis="y", linestyle="--", alpha=alpha)
    plt.tight_layout()
    plt.show()

#############################################################################################
def plot_series_MP(modules, RN, t_span=None, n_steps=None, species_labels=None, title_plot="Metapopulation Dynamics"):
    """
    Graphs the spatiotemporal simulation results of the metapopulation in subplots per species.

    Parameters:
        modules (dict): Simulation results for each patch.
        RN (dict, optional): Reaction network, used to access `SpStr` for species names.
        t_span (tuple, optional): Time range for the simulation (t0, tf). Defaults to None.
        n_steps (int, optional): Number of time steps for the simulation. Defaults to None.
        species_labels (list, optional): List of species names. If None, they are taken from `RN.SpStr` if available.
        title_plot (str, optional): Title for the overall figure. Defaults to "Metapopulation Dynamics".
    
    Raises:
        ValueError: If `modules` is None or empty.
    """

    if modules is None or len(modules) == 0:
        raise ValueError("Simulation results must be provided.")

    if n_steps is None:
        n_steps = 500  # Default number of steps  

    if t_span is None:
        t_span = (0, 100)  # Default simulation time range

    # Use species names from RN.SpStr if not provided
    if species_labels is None and RN is not None:
        species_labels = RN.SpStr
    elif species_labels is None:
        species_labels = [f"Species {i + 1}" for i in range(next(iter(modules.values())).shape[1])]

    t = np.linspace(t_span[0], t_span[1], n_steps)
    n_species = len(species_labels)

    # Create a figure with subplots for each species
    fig, axes = plt.subplots(n_species, 1, figsize=(10, 4 * n_species), sharex=True)

    # Add a general title to the figure
    fig.suptitle(title_plot, fontsize=16, fontweight="bold")

    if n_species == 1:
        axes = [axes]  # Ensure 'axes' is iterable if there is only one species

    for idx, ax in enumerate(axes):
        for patch, values in modules.items():
            ax.plot(t, values[:, idx], label=f"{patch}") 
        ax.set_xlabel("Time")
        ax.set_ylabel(f"{species_labels[idx]}")
        ax.legend(loc="upper right", fontsize="small", ncol=1)
        ax.grid(True)

    plt.tight_layout() 
    plt.show()

# Function to plot the dynamics of all species in 2D
def plot_all_species_2d(results, species_names=None):
    """
    Generates 2D phase diagrams for all combinations of species.

    Parameters:
        results (dict): Output from simulate_odes_metapop_mak with data from each patch.
        species_names (list, optional): Names of species used for labels.

    Returns:
        None
    """
    num_patches = len(results)  # Number of patches
    num_species = results[list(results.keys())[0]].shape[1]  # Number of species

    # Species names
    if species_names is None:
        species_names = [f"Species {i+1}" for i in range(num_species)]

    # Create a figure with subplots for each pair of species
    fig_cols = min(3, num_species)  # Maximum 3 columns
    fig_rows = int(np.ceil(num_species * (num_species - 1) / 2 / fig_cols))  # Required rows

    fig, axes = plt.subplots(fig_rows, fig_cols, figsize=(fig_cols * 6, fig_rows * 6))
    
    # If there is only one row, make sure axes is 1D
    if fig_rows == 1:
        axes = axes.flatten()

    # Iterate over all possible combinations of species
    idx = 0
    for i in range(num_species):
        for j in range(i + 1, num_species):
            x_label = species_names[i]
            y_label = species_names[j]

            # Plot all patch combinations for this species pair
            for patch, data in results.items():
                x = data[:, i]
                y = data[:, j]

                ax = axes[idx] if fig_rows == 1 else axes[idx // fig_cols, idx % fig_cols]
                ax.plot(x, y, label=f"{patch}", alpha=0.7)

            # Adjust titles and labels
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label) 
            ax.legend()

            idx += 1

    # Adjust layout and display the plot
    plt.tight_layout()
    plt.show()

# Function to plot the dynamics of all species in 3D
def plot_all_species_3d(results, species_names=None):
    """
    Generates 3D plots with all possible combinations of three species in the same graph.

    Parameters:
        results (dict): Output from simulate_odes_metapop_mak with data from each patch.
        species_names (list, optional): Names of species used for labels.

    Returns:
        None
    """
    num_species = next(iter(results.values())).shape[1]  # Number of species
    
    if num_species < 3:
        raise ValueError("At least three species are required to generate 3D plots.")
    
    # Generate all possible combinations of three species
    species_combinations = list(itertools.combinations(range(num_species), 3))
    
    # Assign species names if not provided
    if species_names is None:
        species_names = [f"Species {j+1}" for j in range(num_species)]
    
    # Create 3D figures for each species combination
    for comb in species_combinations:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        for patch, data in results.items():
            x, y, z = data[:, comb[0]], data[:, comb[1]], data[:, comb[2]]
            ax.plot(x, y, z, label=f"{patch}", alpha=0.7)
        
        # Labels and title
        ax.set_xlabel(species_names[comb[0]])
        ax.set_ylabel(species_names[comb[1]])
        ax.set_zlabel(species_names[comb[2]])
        ax.set_title(f"3D Dynamics: {species_names[comb[0]]}, {species_names[comb[1]]}, {species_names[comb[2]]}")
        
        ax.legend()
        plt.tight_layout()
        plt.show()

#############################################################################################
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

def plot_series_PDE(simulation_data, species_names, t_span, time_points=None):
    n_species = simulation_data.shape[-1]
    
    if time_points is None:
        time_points = [t_span[0], t_span[1] // 2, t_span[1]]
    
    for s in range(n_species):
        fig, axes = plt.subplots(1, len(time_points), figsize=(15, 5), constrained_layout=True)
        fig.suptitle(f"Concentration of {species_names[s]}")
        
        vmin, vmax = np.min(simulation_data[:, :, :, s]), np.max(simulation_data[:, :, :, s])
        im = None
        
        for idx, t in enumerate(time_points):
            ax = axes[idx]
            im = ax.imshow(simulation_data[t, :, :, s], cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
            ax.set_title(f"Time: {t}")
            ax.set_xticks(np.arange(simulation_data.shape[2]))
            ax.set_yticks(np.arange(simulation_data.shape[1]))
        
        fig.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical', fraction=0.02, pad=0.04)
        plt.show()




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider

def create_animation(simulation_data, species_idx, species_name, interval=100):
    fig, ax = plt.subplots()
    
    # Crear un eje para el deslizador
    ax_slider = plt.axes([0.1, 0.02, 0.8, 0.03], facecolor='lightgoldenrodyellow')
    
    # Inicializar el deslizador con el rango adecuado
    slider = Slider(ax_slider, 'Time', 0, simulation_data.shape[0]-1, valinit=0, valstep=1)

    # Inicializar la imagen con el primer frame
    im = ax.imshow(simulation_data[0, :, :, species_idx], cmap='viridis', origin='lower')
    cbar = fig.colorbar(im, ax=ax)  # Barra de colores

    # Función de actualización para la animación
    def update(frame):
        im.set_array(simulation_data[frame, :, :, species_idx]) # Actualizar la imagen
        ax.set_title(f"{species_name} - Time: {frame}")         # Actualizar el título
        slider.set_val(frame)                # Sincroniza el deslizador con el frame actual
        slider.valtext.set_text(f"{frame}")  # Muestra el número en el deslizador
        return im,

    # Función de actualización para el deslizador
    def slider_update(val):
        frame = int(slider.val)
        update(frame)
        plt.draw()

    # Conectar el deslizador con la función de actualización
    slider.on_changed(slider_update)
    
    # Crear la animación
    ani = animation.FuncAnimation(fig, update, frames=simulation_data.shape[0], interval=interval)

    plt.show()
    return ani

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint

# Suponiendo que la función simulate_pde_rd ya está definida
def create_species_animation(RN, result, interval=50, save=False):
    """
    Genera una animación para cada especie en la simulación y opcionalmente la guarda como archivo .mp4.
    
    Parámetros:
    - RN: Red de reacción con nombres de especies en RN.SpStr
    - result: Array de la simulación con forma (n_steps, grid_x, grid_y, n_species)
    - interval: Tiempo en milisegundos entre cuadros
    - save: Si es True, guarda la animación en formato MP4
    """
    n_steps, grid_x, grid_y, n_species = result.shape
    t_values = np.linspace(0, 100, n_steps)
    
    for i in range(n_species):
        fig, ax = plt.subplots()
        im = ax.imshow(result[0, :, :, i], cmap='viridis', animated=True, vmin=result[:, :, :, i].min(), vmax=result[:, :, :, i].max())
        ax.set_title(f"Evolución de {RN.SpStr[i]}")
        
        def update(frame):
            im.set_array(result[frame, :, :, i])
            return im,
        
        ani = animation.FuncAnimation(fig, update, frames=n_steps, interval=interval, blit=True)
        
        if save:
            ani.save(f"species_{RN.SpStr[i]}.mp4", writer="ffmpeg")
        else:
            plt.show()
        
        plt.close(fig)

# Ejemplo de uso
# RN = definir_red_reaccion()
# result = simulate_pde_rd(RN)
# create_species_animation(RN, result, save=True)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button

def animate_series_PDE(simulation_data, species_names,t_span):
    n_species = simulation_data.shape[-1]
    n_time = t_span[1] 
    
    fig, axes = plt.subplots(1, n_species, figsize=(5 * n_species, 5), constrained_layout=True)
    fig.suptitle("Evolution of Species Concentrations")

    if n_species == 1:
        axes = [axes]
    
    ims = []
    
    # Obtener valores mínimo y máximo para normalizar la escala de color
    vmin, vmax = np.min(simulation_data), np.max(simulation_data)
    
    for s, ax in enumerate(axes):
        im = ax.imshow(simulation_data[0, :, :, s], cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"{species_names[s]}")
        ims.append(im)

        # Asegurar que los ejes solo tengan valores enteros
        ax.set_xticks(np.arange(simulation_data.shape[2]))
        ax.set_yticks(np.arange(simulation_data.shape[1]))

    # Agregar una única barra de colores a la derecha
    cbar = fig.colorbar(ims[0], ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
    
    # Crear slider para el tiempo
    ax_slider = plt.axes([0.2, 0.02, 0.6, 0.03])
    slider = Slider(ax_slider, 'Time', 0, n_time, valinit=0, valstep=1)
    
    def update(frame):
        t = int(slider.val)
        for s in range(n_species):
            ims[s].set_array(simulation_data[t, :, :, s])
        fig.canvas.draw_idle()
    
    slider.on_changed(update)
    plt.show()

# def animate_series_PDE(simulation_data, species_names, t_span,interval):
#     n_species = simulation_data.shape[-1]
#     n_time = t_span[1]
    
#     fig, axes = plt.subplots(1, n_species, figsize=(5 * n_species, 5), constrained_layout=True)
#     fig.suptitle("Evolution of Species Concentrations")

#     if n_species == 1:
#         axes = [axes]
    
#     ims = []
    
#     # Obtener valores mínimo y máximo para normalizar la escala de color
#     vmin, vmax = np.min(simulation_data), np.max(simulation_data)
    
#     for s, ax in enumerate(axes):
#         im = ax.imshow(simulation_data[0, :, :, s], cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
#         ax.set_title(f"{species_names[s]}")
#         ims.append(im)

#         # Asegurar que los ejes solo tengan valores enteros
#         ax.set_xticks(np.arange(simulation_data.shape[2]))
#         ax.set_yticks(np.arange(simulation_data.shape[1]))

#     # Agregar una única barra de colores a la derecha
#     cbar = fig.colorbar(ims[0], ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
    
#     # Crear slider para el tiempo
#     ax_slider = plt.axes([0.1, 0.02, 0.6, 0.03])
#     slider = Slider(ax_slider, 'Time', 0, n_time, valinit=0, valstep=1)
    
#     # Botones Play y Pause
#     ax_play = plt.axes([0.75, 0.02, 0.08, 0.04])
#     ax_pause = plt.axes([0.85, 0.02, 0.08, 0.04])
    
#     btn_play = Button(ax_play, "Play")
#     btn_pause = Button(ax_pause, "Pause")

#     anim_running = [True]  # Variable mutable para controlar la animación

#     def update(frame):
#         """Actualiza la imagen en cada paso de la animación."""
#         t = frame % (n_time + 1)
#         slider.set_val(t)  # Mueve el slider
#         for s in range(n_species):
#             ims[s].set_array(simulation_data[t, :, :, s])
#         return ims

#     ani = animation.FuncAnimation(fig, update, frames=n_time+1, interval=interval, repeat=True)

#     def play(event):
#         """Reanuda la animación."""
#         if not anim_running[0]:
#             ani.event_source.start()
#             anim_running[0] = True

#     def pause(event):
#         """Pausa la animación."""
#         if anim_running[0]:
#             ani.event_source.stop()
#             anim_running[0] = False

#     btn_play.on_clicked(play)
#     btn_pause.on_clicked(pause)

#     plt.show()


######################################################################################
######################################################################################
######################################################################################
# COGNITIVE DOMAIN
######################################################################################
######################################################################################
######################################################################################
# Función para clasificar un vector de proceso v en categorías base y extendidas
def classify_process(v, S, v_prev=None, tol=1e-3):
    """
    Clasifica un vector de proceso v con respecto a la matriz estequiométrica S.

    Categorías base:
    - Stationary Mode: El proceso no produce cambios netos en las concentraciones (Sv=0).
    - Problem: El proceso consume especies pero no produce ninguna (Sv <= 0 y al menos un Sv < 0).
    - Challenge: El proceso consume al menos una especie (al menos un Sv < 0).
    - Cognitive Domain: El proceso mantiene o aumenta todas las especies (Sv >= 0) y utiliza todas las reacciones de la red (v > 0).
    
    Categorías extendidas (requieren un proceso previo 'v_prev'):
    - Counteraction: Un proceso 'v' que, combinado con un 'Challenge' previo 'v_prev', resulta en un no-consumo neto (S(v + v_prev) >= 0).
    - Solution: Un proceso 'v' del espacio 'Challenge' que sirve como solución para un 'Problem' previo 'v_prev' (S(v+v_prev) >= 0).
    - Cognitive Control: Un proceso 'v' del 'Cognitive Domain' que sirve como solución para un 'Problem' previo 'v_prev' (S(v+v_prev) >= 0, v > 0).

    Args:
        v (np.ndarray): El vector de proceso a clasificar. Debe ser un array de NumPy que representa los flujos.
        S (np.ndarray): La matriz estequiométrica, donde las filas son especies y las columnas son reacciones.
        v_prev (np.ndarray, optional): Un vector de proceso previo, utilizado para clasificar
                                     'Counteraction' y 'Cognitive Control'. Por defecto es None.
        tol (float): Tolerancia para comparaciones numéricas con cero para manejar la imprecisión de punto flotante. Por defecto es 1e-8.
        
    Returns:
        list[str]: Una lista ordenada de las categorías a las que pertenece el proceso v.
                   Si no coincide con ninguna categoría específica, se clasifica como "Other".
                   Un proceso puede pertenecer a múltiples categorías.
    
    Raises:
        TypeError: Si 'v', 'S', o 'v_prev' no son arrays de NumPy.
    """
    
    # Validaciones de tipo para los inputs
    if not isinstance(v, np.ndarray) or not isinstance(S, np.ndarray):
        raise TypeError("Los inputs 'v' y 'S' deben ser arrays de NumPy.")
    if v_prev is not None and not isinstance(v_prev, np.ndarray):
        raise TypeError("El input 'v_prev', si se proporciona, debe ser un array de NumPy.")

    # Calcular el cambio neto en las concentraciones de las especies (Sv)
    Sv = S @ v
    classifications = []

    # --- Propiedades fundamentales de Sv y v (calculadas una vez para eficiencia) ---
    is_stationary_mode = np.all((-tol <= Sv) & (Sv <= tol))         # p.all(np.abs(Sv) < tol)
    all_Sv_non_negative = np.all(Sv >= -tol) and np.any(Sv > tol)   # np.all(Sv >= -tol) 
    has_net_consumption = np.any(Sv < -tol)                      # (np.any(Sv < 0) and np.any((0 <= Sv) & (Sv <= tol))) #
    all_Sv_non_positive = np.all(Sv <= tol) and np.any(Sv <= -tol) # np.all(Sv <= tol) 
    all_Sv_negative = np.all(Sv < -tol) and np.all(Sv < -tol)       # Todos los Sv < 0 (nueva condición para Challenge)
    all_reactions_active = np.all(v > tol) 

    # --- Clasificaciones Base (basadas solo en 'v' y 'S') --- 
    if is_stationary_mode:
        classifications = ["Stationary Mode"]
    elif all_Sv_non_negative and all_reactions_active:
        classifications = ["Cognitive Domain"]        
    elif has_net_consumption and all_Sv_non_positive:
        classifications = ["Problem"]
    elif has_net_consumption:
        classifications = ["Challenge"]
    elif all_Sv_non_negative:
        classifications = ["Overproduction Mode"]  
    elif all_Sv_negative:
        classifications = ["Not Feasible"] 
    else:
        classifications = ["Other"]

    # --- Clasificaciones Extendidas (requieren 'v_prev' para el contexto) ---
    if v_prev is not None:
        Sv_prev = S @ v_prev
        is_v_prev_a_challenge = np.any(Sv_prev < -tol) # np.any(Sv < -tol)
        is_v_prev_a_problem = is_v_prev_a_challenge and (np.all(Sv_prev <= 0) and np.any(Sv_prev <= -tol)) #np.any(Sv < -tol) and (np.all(Sv <= 0) and np.any(Sv <= -tol))

        Sv_combined = S @ (v + v_prev)

        # 6. Counteraction (Contramedida) → solo si fue Challenge pero no Problem
        if is_v_prev_a_challenge and not is_v_prev_a_problem and "Challenge":
            if np.all(Sv_combined >= -tol):
                classifications = ["Counteraction"]

        # 7. Solution (Solución extendida)
        elif is_v_prev_a_problem and "Overproduction Mode" in classifications:
            if np.all(Sv_combined >= -tol):
                classifications = ["Solution"]

        # 8. Cognitive Control (Control Cognitivo)
        elif is_v_prev_a_problem and "Cognitive Domain" in classifications:
            if np.any(Sv_combined >= -tol):
                classifications = ["Cognitive Control"]

    # Devolver una lista de clasificaciones únicas y ordenadas para una salida consistente.
    return sorted(list(set(classifications)))

#####################################################################################
# Función para graficar el histograma de tipos de proceso y guardar en Excel
def plot_process_types_histogram(flux_vector, S, 
                                xlabel="Tipo de Proceso", ylabel="Frecuencia",
                                title="Histograma de Tipos de Proceso", 
                                excel_filename="classified_processes.xlsx",
                                filename="histograma.png", 
                                save_figure=True, 
                                ax=None, show_fig=False):
    """
    Clasifica, grafica histograma, y guarda resultados en Excel.
    """
    # if not isinstance(flux_vector, pd.DataFrame):
    #     raise TypeError("flux_vector debe ser un DataFrame de pandas.")
    if isinstance(flux_vector, np.ndarray):
        flux_vector = pd.DataFrame(flux_vector, columns=[f"v{i+1}" for i in range(flux_vector.shape[1])])
        flux_vector.insert(0, "Time", range(len(flux_vector)))
    elif not isinstance(flux_vector, pd.DataFrame):
        raise TypeError("flux_vector debe ser un DataFrame o un ndarray de NumPy.")    
    if not isinstance(S, np.ndarray):
        raise TypeError("S debe ser un array de NumPy.")    

    # --- Clasificar ---
    flux_values = flux_vector.iloc[:, 1:]  # quitar columna Time
    process_types = []
    for _, row in flux_values.iterrows():
        v = row.to_numpy()
        cat = classify_process(v, S)  # ahora solo se pasa v y S
        # Convertir la lista de clasificaciones a una cadena para hacerla hashable
        cat_str = ",".join(cat) if cat else "None" # ",".join(sorted(cat)) if cat else "None"
        process_types.append(cat_str)

    # Calcular S*v
    Sv_matrix = flux_values.apply(lambda v: S @ v.to_numpy(), axis=1)
    Sv_expanded = pd.DataFrame(Sv_matrix.tolist(), 
                               columns=[f"S*v_{i+1}" for i in range(S.shape[0])],
                               index=flux_vector.index)

    classified_df = pd.concat([flux_vector, Sv_expanded], axis=1)
    classified_df["Process_Type"] = process_types
    # print("classified_df=\n",classified_df)

    # --- Contar ---
    process_counts = Counter(process_types)
    category_order = ["Stationary Mode", "Cognitive Domain", "Problem", "Challenge","Overproduction Mode", "Not Feasible", "Other", "None"]
    labels = [cat for cat in category_order if cat in process_counts]
    counts = [process_counts[cat] for cat in labels]

    # Colores 
    color_map = { 
        "Stationary Mode": "cyan",
        "Cognitive Domain": "green",
        "Problem": "red",
        "Challenge": "orange",
        "Overproduction Mode": "blue",
        "Not Feasible": "yellow",
        "Other": "grey",
        "None": "black"        
    }

    colors = [color_map.get(label, "grey") for label in labels]

    # --- Graficar ---
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.get_figure()

    bars = ax.bar(labels, counts, color=colors)
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval),
                va='bottom', ha='center', fontweight='bold') 

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # --- Guardar en Excel ---
    out_dir = "./visualizations/process_classification"
    os.makedirs(out_dir, exist_ok=True)
    filepath_excel = os.path.join(out_dir, excel_filename)

    with pd.ExcelWriter(filepath_excel, engine="openpyxl") as writer:
        classified_df.to_excel(writer, sheet_name="All_Processes", index=False)
        for category, group in classified_df.groupby("Process_Type", sort=False):
            group.to_excel(writer, sheet_name=category[:31], index=False)
    
    print(f"Procesos clasificados guardados en: {filepath_excel}")

    # --- Guardar figura ---
    if save_figure:
        filepath_fig = os.path.join(out_dir, filename)
        fig.savefig(filepath_fig, dpi=300, bbox_inches="tight")
        print(f"Histograma guardado en: {filepath_fig}")

    if show_fig:
        plt.show()

    return fig, ax, classified_df

#####################################################################################
# Función para graficar el cono y la región factible en 3D
def plot_cone_and_region(S, 
                         grid_max=None, grid_res=5,
                         axis_names=None, show=True,
                         extra_vector=None, 
                         extra_vector_labels=None,
                         extra_vector_colors=None):   
    """
    Genera proyecciones 3D del cono definido por vectores en null_vectors
    y la región factible Sv > 0 para una matriz estequiométrica S.
    Los puntos factibles se clasifican con classify_process() y se pintan
    con diferentes colores según su categoría.
    Además guarda en un archivo Excel los puntos factibles (v), S*v y
    Process_Type en una hoja principal, y en hojas separadas por categoría.
    """
    n = S.shape[1]
    m = S.shape[0]
    if axis_names is None:
        axis_names = [f"v{i+1}" for i in range(n)]
    
    if extra_vector_labels is None:
        extra_vector_labels = [f"Extra vector {i+1}" for i in range(len(extra_vector or []))]
    
    if extra_vector_colors is None:
        extra_vector_colors = ['orange', 'cyan', 'purple', 'brown', 'yellow', 'pink'] + \
                              ['gray'] * max(0, len(extra_vector or []) - 6)
    
    out_dir = "./visualizations/cone_projections_3D"
    os.makedirs(out_dir, exist_ok=True)

    # ---- 1. Construcción del cono (Sv=0) ----
    S_sym = sp.Matrix(S) 
    null_basis = S_sym.nullspace() 
    null_vectors = [np.array(v, dtype=float).flatten() for v in null_basis]
    print("Vectores del espacio nulo =", null_vectors)

    base_max = np.max([np.max(np.abs(v)) for v in null_vectors]) if null_vectors else 1.0

    if extra_vector is not None:
        if isinstance(extra_vector, np.ndarray) and extra_vector.ndim == 1:
            extra_vector = [extra_vector]  
        extra_max = np.max(np.abs(np.vstack(extra_vector))) if extra_vector else 0
    else:
        extra_max = 0

    if grid_max is None:
        grid_max = max(base_max, extra_max)

    # ---- 2. Construcción de la región factible Sv>0 ----
    grid = np.linspace(0, grid_max, grid_res)
    V = np.array(np.meshgrid(*([grid] * n))).T.reshape(-1, n)
    mask = np.all(S @ V.T >= -1e-10, axis=0) 
    points_pos = V[mask]

    # ---- Guardar puntos factibles v, S*v y Process_Type en Excel ----
    if points_pos.shape[0] > 0:
        Sv = (S @ points_pos.T).T  
        df_points = pd.DataFrame(points_pos, columns=axis_names)
        df_Sv = pd.DataFrame(Sv, columns=[f"S*v=x{i+1}" for i in range(m)])

        # Clasificación de procesos
        classifications = []
        for v in points_pos:
            cat = classify_process(v, S)
            cat_str = ",".join(cat) if cat else "None"
            classifications.append(cat_str)
        df_class = pd.DataFrame({"Process_Type": classifications})

        # Combinar todo
        df_all = pd.concat([df_points, df_Sv, df_class], axis=1)

        # Guardar en un solo archivo con varias hojas
        excel_path = os.path.join(out_dir, "classified_processes.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            # Hoja principal con todo
            df_all.to_excel(writer, sheet_name="All_Processes", index=False)

            # Hojas separadas por categoría
            for category, group in df_all.groupby("Process_Type"):
                safe_name = str(category)[:31]  # Excel limita nombres a 31 caracteres
                group.to_excel(writer, sheet_name=safe_name, index=False)

        print(f"Vectores de procesos guardados en: {excel_path}")
    else:
        print("No se encontraron puntos factibles Sv>0.")

    # ---- 3. Proyecciones 3D ----
    saved_files = []
    for (i, j, k) in itertools.combinations(range(n), 3):
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111, projection='3d')

        # Cono (Sv=0)
        if len(null_vectors) == 0:
            ax.scatter([0], [0], [0], color='lightblue', s=50, label="Sv=0 (origen)")
        elif len(null_vectors) == 1:
            v1 = null_vectors[0]
            ax.quiver(0, 0, 0, *(v1[[i,j,k]]), color='black', linewidth=2, 
                      arrow_length_ratio=0.1, label="Stationary mode 1")
        else:
            v1, v2 = null_vectors[:2] 
            tri = np.array([[0,0,0], v1[[i,j,k]], v2[[i,j,k]]])
            poly = Poly3DCollection([tri], alpha=0.4, facecolor='lightblue', 
                                   edgecolor='lightblue', label="Sv=0")
            ax.add_collection3d(poly)
            colors = ['black', 'lightblue', 'red', 'orange', 'purple', 'brown'] + \
                     ['purple'] * (len(null_vectors) - 2)
            for idx, v in enumerate(null_vectors):
                ax.quiver(0, 0, 0, *(v[[i,j,k]]), color=colors[idx], linewidth=2, 
                          arrow_length_ratio=0.1, label=f"Stationary mode {idx+1}")

        # Región factible Sv>0 clasificada
        if points_pos.shape[0] > 0:
            proj_pos = points_pos[:, [i, j, k]]
            color_map = {
                "Stationary Mode": "cyan",
                "Cognitive Domain": "green",
                "Problem": "red",
                "Challenge": "orange",
                "Overproduction Mode": "blue",
                "Not Feasible": "yellow",
                "Other": "grey",
                "None": "black"
            }
            categories = [c.split(",")[0] if c else "None" for c in classifications]

            for cat in sorted(set(categories)):
                mask = [c == cat for c in categories]
                ax.scatter(proj_pos[mask, 0],
                           proj_pos[mask, 1],
                           proj_pos[mask, 2],
                           alpha=0.6, s=20,
                           c=color_map.get(cat, "grey"),
                           label=cat)

        if extra_vector is not None:
            for k_idx, vec in enumerate(extra_vector):
                extra_proj = vec[[i, j, k]]
                ax.quiver(0, 0, 0, *extra_proj,
                          color=extra_vector_colors[k_idx % len(extra_vector_colors)],
                          linewidth=2, arrow_length_ratio=0.1,
                          label=extra_vector_labels[k_idx] if k_idx < len(extra_vector_labels) else f"Extra vector {k_idx+1}")

        ax.set_xlim([-grid_max*0.1, grid_max*1.1])
        ax.set_ylim([-grid_max*0.1, grid_max*1.1])
        ax.set_zlim([-grid_max*0.1, grid_max*1.1])

        ax.set_xlabel(axis_names[i])
        ax.set_ylabel(axis_names[j])
        ax.set_zlabel(axis_names[k])
        ax.set_title(f"Proyección 3D ({axis_names[i]}, {axis_names[j]}, {axis_names[k]})")
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

        filename = os.path.join(out_dir, f"Cone_and_region_{i+1}_{j+1}_{k+1}.png")
        fig.savefig(filename, dpi=150, bbox_inches='tight')
        saved_files.append(filename)
        if show:
            plt.show()
        else:
            plt.close(fig)

    print(f"Se guardaron {len(saved_files)} imágenes en: {out_dir}")
    return saved_files, points_pos

#####################################################################################
# Funciones para graficar series de tiempo con intervalos de Cognitive Control
def plot_series_with_domain_intervals(time_series, flux_vector, S,
                                      title="Serie de Tiempo de Concentraciones",
                                      save_figure=False,
                                      ax=None,
                                      show_fig=False):
    fig, ax = plot_series_ode(time_series, title=title, save_figure=save_figure, ax=ax)

    times = time_series["Time"].to_numpy()
    flux_values = flux_vector.iloc[:, 1:]
    process_types = [classify_process(v.to_numpy(), S) for _, v in flux_values.iterrows()]

    def get_intervals(mask, times):
        intervals = []
        in_interval = False
        start = None
        for t, flag in zip(times, mask):
            if flag and not in_interval:
                in_interval = True
                start = t
            elif not flag and in_interval:
                in_interval = False
                intervals.append((start, t))
        if in_interval:
            intervals.append((start, times[-1]))
        return intervals

    is_cd = np.array(["Cognitive Domain" in cat for cat in process_types])
    is_sm = np.array(["Stationary Mode" in cat for cat in process_types])

    cd_intervals = get_intervals(is_cd, times)
    sm_intervals = get_intervals(is_sm, times)

    for (t_start, t_end) in cd_intervals:
        ax.axvspan(t_start, t_end, color="darkgreen", alpha=0.2)
        ax.axvline(x=t_start, color="darkgreen", linestyle="--", alpha=0.8)
        ax.axvline(x=t_end, color="darkgreen", linestyle="--", alpha=0.8)

    for (t_start, t_end) in sm_intervals:
        ax.axvspan(t_start, t_end, color="cyan", alpha=0.15)
        ax.axvline(x=t_start, color="cyan", linestyle="--", alpha=0.7)
        ax.axvline(x=t_end, color="cyan", linestyle="--", alpha=0.7)

    ax.plot([], [], color="darkgreen", linestyle="--", label="Cognitive Control")
    ax.plot([], [], color="cyan", linestyle="--", label="Stationary Mode")
    ax.legend(loc="upper right")

    if show_fig:
        plt.show()

    return fig, ax

# Función para graficar flujos con intervalos de Cognitive Domain
def plot_flux_with_domain_intervals(flux_vector, S,
                                    title="Serie de Tiempo de Flujos",
                                    save_figure=False,
                                    ax=None,
                                    show_fig=False):
    """
    Extiende plot_series_ode aplicada al flux_vector para resaltar los intervalos
    donde ocurren 'Cognitive Domain' (verde) y 'Stationary Mode' (rojo).
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Graficar la serie original
    fig, ax = plot_series_ode(flux_vector, title=title, save_figure=save_figure, ax=ax)
    print("\nflux_vector_shape =", flux_vector.shape)
    
    # Clasificar procesos
    flux_values = flux_vector.iloc[:, 1:]  # quitar columna Time
    process_types = []
    for _, row in flux_values.iterrows():
        v = row.to_numpy()
        cat = classify_process(v, S)
        process_types.append(cat)
    
    # Extraer tiempos
    times = flux_vector["Time"].to_numpy()

    # Función auxiliar para obtener intervalos consecutivos True
    def get_intervals(mask, times):
        intervals = []
        in_interval = False
        start = None
        for t, flag in zip(times, mask):
            if flag and not in_interval:
                in_interval = True
                start = t
            elif not flag and in_interval:
                in_interval = False
                intervals.append((start, t))
        if in_interval:
            intervals.append((start, times[-1]))
        return intervals

    # === 1. Intervalos de Cognitive Domain (verde) ===
    is_cd = np.array(["Cognitive Domain" in cat for cat in process_types])
    cd_intervals = get_intervals(is_cd, times)

    # === 2. Intervalos de Stationary Mode (rojo) ===
    is_sm = np.array(["Stationary Mode" in cat for cat in process_types])
    sm_intervals = get_intervals(is_sm, times)

    # --- Graficar intervalos ---
    for (t_start, t_end) in cd_intervals:
        ax.axvspan(t_start, t_end, color="darkgreen", alpha=0.2)
        ax.axvline(x=t_start, color="darkgreen", linestyle="--", alpha=0.8)
        ax.axvline(x=t_end, color="darkgreen", linestyle="--", alpha=0.8)

    for (t_start, t_end) in sm_intervals:
        ax.axvspan(t_start, t_end, color="cyan", alpha=0.15)
        ax.axvline(x=t_start, color="cyan", linestyle="--", alpha=0.7)
        ax.axvline(x=t_end, color="cyan", linestyle="--", alpha=0.7)

    # --- Estadísticas de duración ---
    def print_stats(name, intervals):
        if len(intervals) == 0:
            print(f"No se detectaron intervalos de {name}.")
            return
        durations = [t2 - t1 for (t1, t2) in intervals]
        print(f"\n{name}: {len(intervals)} intervalos")
        print(f"  T_min = {min(durations):.4f}")
        print(f"  T_max = {max(durations):.4f}")
        print(f"  T_avg = {np.mean(durations):.4f}")

    print_stats("Cognitive Domain", cd_intervals)
    print_stats("Stationary Mode", sm_intervals)

    # Leyenda
    ax.plot([], [], color="darkgreen", linestyle="--", label="Cognitive Control")
    ax.plot([], [], color="cyan", linestyle="--", label="Stationary Mode")
    ax.legend(loc="upper right")

    if show_fig:
        plt.show()

    return fig, ax


#############################################################################################################
# Función para sumar vectores de flux1 y flux2, clasificar y graficar histograma
def histogram_flux_sum(S, flux1, flux2, 
                       title="Histograma de Tipos de Proceso (Suma Flux1 + Flux2)", 
                       filename="histograma_sum.png", csv_filename="sum_flux_data.csv",
                       excel_filename="sum_flux_data.xlsx", save_figure=True, show_fig=True, 
                       max_combinations=None):
    """
    Suma cada vector de flux1 con cada vector de flux2, clasifica los vectores resultantes,
    genera un histograma de tipos de proceso y guarda los resultados en archivos CSV y Excel.
    - CSV: Incluye v_f1, v_f2, v_combined, S*v y Process_Type en el orden automático.
    - Excel: Incluye solo v_combined, S*v y Process_Type, con hojas por Process_Type.

    Parameters:
    - S: Matriz de transformación (numpy array).
    - flux1: DataFrame o NumPy array con columnas Flux_r* y opcionalmente Time (por ejemplo, vectores de desafío).
    - flux2: DataFrame o NumPy array con columnas Flux_r* y opcionalmente Time (por ejemplo, vectores de control cognitivo).
    - title: Título del histograma.
    - filename: Nombre del archivo para guardar el histograma.
    - csv_filename: Nombre del archivo CSV para guardar todos los datos.
    - excel_filename: Nombre del archivo Excel para guardar v_combined, S*v y Process_Type.
    - save_figure: Booleano para guardar la figura.
    - show_fig: Booleano para mostrar la figura.
    - max_combinations: Límite opcional para el número de combinaciones a generar (default: None).

    Returns:
    - fig, ax, combined_df: Objetos de la figura, ejes de Matplotlib y DataFrame con los vectores sumados.
    """
    # Convertir flux1 y flux2 a DataFrame si son NumPy arrays
    if isinstance(flux1, np.ndarray):
        flux1 = pd.DataFrame(flux1, columns=[f'Flux_r{i+1}' for i in range(flux1.shape[1])])
    if isinstance(flux2, np.ndarray):
        flux2 = pd.DataFrame(flux2, columns=[f'Flux_r{i+1}' for i in range(flux2.shape[1])])

    # Extraer las columnas de flujos
    flux_columns = [col for col in flux1.columns if col.startswith('Flux_r')]
    if not all(col in flux2.columns for col in flux_columns):
        raise ValueError("flux1 y flux2 deben tener las mismas columnas de flujos (Flux_r*).")

    print(f"Se encontraron {len(flux1)} vectores en flux1.")
    print(f"Se encontraron {len(flux2)} vectores en flux2.")

    # Verificar que ambos DataFrames no estén vacíos
    if flux1.empty or flux2.empty:
        raise ValueError("Uno o ambos DataFrames (flux1 o flux2) están vacíos.")

    # Verificar si las columnas 'Time' están presentes
    has_time_f1 = 'Time' in flux1.columns
    has_time_f2 = 'Time' in flux2.columns
    print(f"flux1 tiene columna 'Time': {has_time_f1}")
    print(f"flux2 tiene columna 'Time': {has_time_f2}")

    # Generar combinaciones sumadas de vectores
    flux_data = []
    combination_count = 0
    for idx_f2, row_f2 in flux2.iterrows():
        time_f2 = row_f2['Time'] if has_time_f2 else idx_f2
        v_f2 = row_f2[flux_columns].to_numpy()
        for idx_f1, row_f1 in flux1.iterrows():
            if max_combinations is not None and combination_count >= max_combinations:
                break
            time_f1 = row_f1['Time'] if has_time_f1 else idx_f1
            v_f1 = row_f1[flux_columns].to_numpy()
            # Sumar los vectores
            v_combined = v_f2 + v_f1
            # Usar el tiempo de flux2 o promedio si ambos tienen 'Time'
            time_value = time_f2 if not has_time_f1 else (time_f1 + time_f2) / 2 if has_time_f2 else time_f1
            # Crear diccionario para el vector combinado
            row_data = {}
            # Agregar v_combined
            for col, val in zip(flux_columns, v_combined):
                row_data[col] = val
            # Agregar v_f1
            for i, val in enumerate(v_f1):
                row_data[f'Flux1_r{i+1}'] = val
            # Agregar v_f2
            for i, val in enumerate(v_f2):
                row_data[f'Flux2_r{i+1}'] = val
            # Agregar Time si está presente
            if has_time_f1 or has_time_f2:
                row_data['Time'] = time_value
            flux_data.append(row_data)
            combination_count += 1
        if max_combinations is not None and combination_count >= max_combinations:
            break

    # Crear DataFrame con los vectores combinados
    flux_vector = pd.DataFrame(flux_data)
    # flux_vector_df = flux_vector
    # print("flux_vector=\n", flux_vector_df)
    print(f"Total de vectores combinados generados: {len(flux_vector)}")

    # Clasificar los procesos en el flux_vector combinado
    flux_values = flux_vector[flux_columns]  # Excluir columna Time si existe
    process_types = []
    Sv_values = []
    for _, row in flux_values.iterrows():
        v = row.to_numpy()
        # Usar la lógica de classify_process
        Sv = S @ v
        Sv_values.append(Sv) 
        if np.all((0 < Sv) & (Sv <= 1e-8)):
            category = "Stationary mode"
        elif np.all(Sv > 0) and np.any(Sv > 1e-8):
            category = "Cognitive Control"
        elif np.any(Sv < 0) and np.any(Sv > 0):
            category = "Challenge"
        elif np.any(Sv < 0):
            category = "Problem"
        else:
            category = "Other"
        process_types.append(category)

    # Calcular S*v
    Sv_expanded = pd.DataFrame(Sv_values, 
                               columns=[f"S*v_{i+1}" for i in range(S.shape[0])],
                               index=flux_vector.index)

    # Crear DataFrame combinado
    combined_df = pd.concat([flux_vector, Sv_expanded], axis=1)
    combined_df["Process_Type"] = process_types

    # Determinar automáticamente el orden de las columnas para CSV
    csv_columns = []
    if has_time_f1 or has_time_f2:
        csv_columns.append('Time')
    csv_columns.extend([col for col in combined_df.columns if col.startswith('Flux1_r')])
    csv_columns.extend([col for col in combined_df.columns if col.startswith('Flux2_r')])
    csv_columns.extend([col for col in combined_df.columns if col.startswith('Flux_r')])
    csv_columns.extend([col for col in combined_df.columns if col.startswith('S*v_')])
    csv_columns.append('Process_Type')
    csv_available_columns = [col for col in csv_columns if col in combined_df.columns]
    csv_df = combined_df.reindex(columns=csv_available_columns)
    # print("CSV columns=\n", csv_df.columns)

    # Determinar columnas para Excel (v_combined, S*v, Process_Type)
    excel_columns = []
    if has_time_f1 or has_time_f2:
        excel_columns.append('Time')
    excel_columns.extend([col for col in combined_df.columns if col.startswith('Flux_r')])
    excel_columns.extend([col for col in combined_df.columns if col.startswith('S*v_')])
    excel_columns.append('Process_Type')
    excel_available_columns = [col for col in excel_columns if col in combined_df.columns]
    excel_df = combined_df.reindex(columns=excel_available_columns)
    # print("Excel columns=\n", excel_df.columns)

    # Contar frecuencias para el histograma
    process_counts = Counter(process_types)
    category_order = ["Stationary mode", "Cognitive Control", "Problem", "Challenge", "Other"]
    labels = [cat for cat in category_order if cat in process_counts]
    counts = [process_counts[cat] for cat in labels]

    # Colores para el histograma 
    color_map = { 
        "Stationary mode": "cyan",
        "Cognitive Control": "green",
        "Problem": "red",
        "Challenge": "orange",
        "Overproduction Mode": "blue",
        "Not Feasible": "yellow",
        "Other": "grey",
        "None": "black"        
    }
    colors = [color_map.get(label, "grey") for label in labels]

    # Graficar
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(labels, counts, color=colors)
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval),
                va='bottom', ha='center', fontweight='bold')

    ax.set_xlabel("Tipo de Proceso", fontsize=12)
    ax.set_ylabel("Frecuencia", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Guardar en CSV
    out_dir = "./visualizations/process_classification/data_xslx"
    os.makedirs(out_dir, exist_ok=True)

    # Definir la subcarpeta de salida
    out_dir_sub = "./visualizations/process_classification/data_csv"
    os.makedirs(out_dir_sub, exist_ok=True)  # Crear la subcarpeta si no existe

    try:
        # Guardar el archivo CSV principal en la subcarpeta
        filepath_csv = os.path.join(out_dir_sub, csv_filename)
        csv_df.to_csv(filepath_csv, index=False)
        print(f"Procesos clasificados guardados en CSV: {filepath_csv}")

        # Guardar archivos CSV separados por Process_Type en la subcarpeta
        for category, group in csv_df.groupby("Process_Type"):
            # Asegurar que el nombre del archivo sea válido (reemplazar espacios y limitar longitud)
            safe_category = category[:31].replace(' ', '_').replace('/', '_').replace('\\', '_')
            category_filepath = os.path.join(out_dir_sub, f"{csv_filename[:-4]}_{safe_category}.csv")
            group.to_csv(category_filepath, index=False)
            print(f"Procesos de tipo '{category}' guardados en CSV: {category_filepath}")

    except PermissionError as e:
        warnings.warn(f"Error de permisos al guardar el archivo CSV: {e}")
        print("No se pudo guardar el archivo CSV. Por favor, verifica los permisos de escritura en la carpeta.")
    except OSError as e:
        warnings.warn(f"Error del sistema al guardar el archivo CSV: {e}")
        print("No se pudo guardar el archivo CSV. Por favor, verifica el espacio en disco o la validez del nombre del archivo.")
    except Exception as e:
        warnings.warn(f"Error inesperado al guardar el archivo CSV: {e}")
        print("No se pudo guardar el archivo CSV. Ocurrió un error inesperado.")

    # Guardar en Excel
    filepath_excel = os.path.join(out_dir, excel_filename)
    # try:
    with pd.ExcelWriter(filepath_excel, engine="openpyxl") as writer:
        excel_df.to_excel(writer, sheet_name="All_Processes", index=False)
        for category, group in excel_df.groupby("Process_Type"):
            group.to_excel(writer, sheet_name=category[:31], index=False)
    print(f"Procesos clasificados guardados en Excel: {filepath_excel}") 

    # Guardar figura
    out_dir_hist = "./visualizations/process_classification/histogram_flux_sum"
    os.makedirs(out_dir_hist, exist_ok=True)
    if save_figure:
        filepath_fig = os.path.join(out_dir_hist, filename)
        fig.savefig(filepath_fig, dpi=300, bbox_inches="tight")
        print(f"Histograma guardado en: {filepath_fig}")

    if show_fig:
        plt.show()

    return fig, ax, combined_df

#################################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from collections import Counter

def analyze_time_scales_all(flux_vector, S, max_scale, tol=1e-8, show_fig=True):
    """
    Analiza la serie temporal de flujos en diferentes escalas de tiempo para encontrar
    la frecuencia de todos los tipos de procesos.

    Args:
        flux_vector (pd.DataFrame): DataFrame de la simulación. La primera columna
                                    debe ser el tiempo y el resto los flujos.
        S (np.ndarray): La matriz estequiométrica.
        max_scale (int): El tamaño máximo de la ventana de tiempo a analizar. 
        tol (float): Tolerancia numérica para clasificar procesos.
        show_fig (bool): Si True, muestra los gráficos.

    Returns:
        dict: Diccionario anidado con las escalas de tiempo como claves y un
              Counter de categorías de procesos como valores.
    """
    if not isinstance(flux_vector, pd.DataFrame):
        raise TypeError("flux_vector debe ser un DataFrame de pandas.")
    if max_scale < 1 or not isinstance(max_scale, int):
        raise ValueError("max_scale debe ser un entero mayor o igual a 1.")

    print(f"Analizando {max_scale} escalas de tiempo (todas las categorías). Esto puede tardar un momento...")

    flux_values = flux_vector.iloc[:, 1:].to_numpy()  # Excluir columna Time
    results = {}

    # Iterar sobre cada escala de tiempo
    for scale in range(1, max_scale + 1):
        categories_counter = Counter()
        num_windows = len(flux_values) - scale + 1
        
        for i in range(num_windows):
            window = flux_values[i : i + scale]
            v_combined = np.sum(window, axis=0)
            
            # Clasificar el proceso combinado
            classifications = classify_process(v_combined, S, tol=tol)
            
            # Actualizar el contador con todas las categorías asignadas
            categories_counter.update(classifications)

        results[scale] = categories_counter
        print(f"Escala {scale}: {dict(categories_counter)}")

    # --- Visualización de Resultados ---
    if show_fig:
        # Extraer todas las categorías posibles
        all_categories = sorted(set(cat for c in results.values() for cat in c.keys()))

        plt.style.use('seaborn-v0_8-whitegrid')
        plt.figure(figsize=(12, 7))

        for category in all_categories:
            counts = [results[scale].get(category, 0) for scale in results.keys()]
            plt.plot(results.keys(), counts, marker='o', linestyle='-', label=category)

        plt.title('Análisis de la Escala de Tiempo (todas las categorías)', fontsize=16, fontweight='bold')
        plt.xlabel('Escala de Tiempo (Tamaño de la Ventana)', fontsize=12)
        plt.ylabel('Frecuencia de Procesos', fontsize=12)
        plt.xticks(list(results.keys()))
        plt.legend()

        # Guardar la figura
        out_dir = "./visualizations/process_classification"
        filename = "time_scale_analysis_all.png"
        os.makedirs(out_dir, exist_ok=True)
        save_path = os.path.join(out_dir, filename)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figura guardada en: {save_path}")

        plt.show()

    return results

def analyze_time_scales_selected(flux_vector, S, max_scale, categories_to_plot=None, tol=1e-8, show_fig=True):
    """
    Analiza la serie temporal de flujos en diferentes escalas de tiempo
    y permite graficar solo las categorías seleccionadas.

    Args:
        flux_vector (pd.DataFrame): DataFrame de la simulación. La primera columna debe ser el tiempo.
        S (np.ndarray): La matriz estequiométrica.
        max_scale (int): El tamaño máximo de la ventana de tiempo a analizar.
        categories_to_plot (list, optional): Lista de categorías a graficar. Si None, grafica todas.
        tol (float): Tolerancia para clasificar procesos.
        show_fig (bool): Si True, muestra los gráficos.

    Returns:
        dict: Diccionario con las frecuencias de procesos por escala.
    """
    flux_values = flux_vector.iloc[:, 1:].to_numpy()
    results = {}

    for scale in range(1, max_scale + 1):
        categories_counter = Counter()
        num_windows = len(flux_values) - scale + 1
        for i in range(num_windows):
            window = flux_values[i : i + scale]
            v_combined = np.sum(window, axis=0)
            classifications = classify_process(v_combined, S, tol=tol)
            categories_counter.update(classifications)
        results[scale] = categories_counter

    if show_fig:
        # Determinar categorías a graficar
        if categories_to_plot is None:
            all_categories = sorted(set(cat for c in results.values() for cat in c.keys()))
        else:
            all_categories = categories_to_plot

        plt.style.use('seaborn-v0_8-whitegrid')
        plt.figure(figsize=(12, 7))

        for category in all_categories:
            counts = [results[scale].get(category, 0) for scale in results.keys()]
            plt.plot(results.keys(), counts, marker='o', linestyle='-', label=category)

        plt.title('Análisis de la Escala de Tiempo (categorías seleccionadas)', fontsize=16, fontweight='bold')
        plt.xlabel('Escala de Tiempo (Tamaño de la Ventana)', fontsize=12)
        plt.ylabel('Frecuencia de Procesos', fontsize=12)
        plt.xticks(list(results.keys()))
        plt.legend()

        out_dir = "./visualizations/process_classification"
        os.makedirs(out_dir, exist_ok=True)
        save_path = os.path.join(out_dir, "time_scale_analysis_selected.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figura guardada en: {save_path}")

        plt.show()

    return results


#################################################################################################
    
import os
import matplotlib.pyplot as plt
from collections import Counter

def analyze_scales(flux_vector, S, max_group_size=5, tol=1e-8, show_fig=True, save_fig=True):
    """
    Analiza la serie temporal de procesos en diferentes escalas de agrupación
    y grafica la frecuencia del 'Cognitive Domain' según la escala.
    """
    if not isinstance(flux_vector, pd.DataFrame):
        raise TypeError("flux_vector debe ser un DataFrame de pandas.")
    if max_group_size < 1 or not isinstance(max_group_size, int):
        raise ValueError("max_group_size debe ser un entero mayor o igual a 1.")

    print(f"Analizando {max_group_size} escalas de tiempo. Esto puede tardar un momento...")

    # Convertir a numpy y trasponer: (n_reacciones x n_tiempo)
    flux_vector = flux_vector.iloc[:, 1:].to_numpy().T  

    results = {}
    for group_size in range(1, max_group_size+1):
        classifications = []

        for i in range(0, flux_vector.shape[1] - group_size + 1, group_size):
            v_combined = np.sum(flux_vector[:, i:i+group_size], axis=1)  # (n_reacciones,)
            classes = classify_process(v_combined, S, tol=tol)
            classifications.extend(classes)

        # Guardar conteo por categoría
        results[group_size] = dict(Counter(classifications))

    # ================================
    # VISUALIZACIÓN SOLO PARA "Cognitive Domain"
    # ================================
    if show_fig and results:
        scales = list(results.keys()) # Escalas de tiempo
        counts = [results[s].get("Cognitive Domain", 0) for s in scales] # Conteos

        # Encontrar escala óptima
        optimal_scale = scales[int(np.argmax(counts))]
        max_count = max(counts)

        plt.style.use('seaborn-v0_8-whitegrid')
        plt.figure(figsize=(12, 7))
        plt.plot(scales, counts, marker='o', linestyle='-', color='b', label='Conteo de Dominio Cognitivo')

        # Resaltar el óptimo
        plt.axvline(x=optimal_scale, color='r', linestyle='--', 
                    label=f'Escala Óptima = {optimal_scale} (Conteo: {max_count})')
        plt.scatter(optimal_scale, max_count, color='red', s=100, zorder=5)

        plt.title('Análisis de la Escala de Tiempo de Estabilidad', fontsize=16, fontweight='bold')
        plt.xlabel('Escala de Tiempo (Tamaño de la Ventana)', fontsize=12)
        plt.ylabel('Frecuencia de Procesos de "Dominio Cognitivo"', fontsize=12)
        plt.xticks(scales)
        plt.legend()

        # Guardar automáticamente
        if save_fig:
            out_dir = "./visualizations/process_classification"
            filename = "time_scale_analysis.png"
            os.makedirs(out_dir, exist_ok=True)
            save_path = os.path.join(out_dir, filename)
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            print(f"Figura guardada en: {save_path}")

        plt.show()

    return results


#################################################################################################
def analyze_time_scales(flux_vector, S, max_scale):
    """
    Analiza la serie temporal de flujos en diferentes escalas de tiempo para encontrar
    la escala óptima de automantenimiento (Dominio Cognitivo).

    Args:
        flux_vector (pd.DataFrame): DataFrame de la simulación. La primera columna
                                    debe ser el tiempo y el resto los flujos.
        S (np.ndarray): La matriz estequiométrica.
        max_scale (int): El tamaño máximo de la ventana de tiempo a analizar. 

    Returns:
        dict: Un diccionario con las escalas de tiempo como claves y el conteo de
              procesos de "Dominio Cognitivo" como valores.
    """
    if not isinstance(flux_vector, pd.DataFrame):
        raise TypeError("flux_vector debe ser un DataFrame de pandas.")
    if max_scale < 1 or not isinstance(max_scale, int):
        raise ValueError("max_scale debe ser un entero mayor o igual a 1.")

    print(f"Analizando {max_scale} escalas de tiempo. Esto puede tardar un momento...")

    # Extraer solo los valores de flujo como un array de NumPy para eficiencia
    flux_values = flux_vector.iloc[:, 1:].to_numpy() # Excluir columna Time
    
    # Diccionario para almacenar los resultados
    results = {}

    # Iterar sobre cada escala de tiempo (tamaño de la ventana)
    for scale in range(1, max_scale + 1):
        cognitive_domain_count = 0
        
        # Deslizar la ventana a través de la serie temporal
        # El bucle se detiene antes para asegurar que la ventana no exceda los límites
        num_windows = len(flux_values) - scale + 1
        for i in range(num_windows):
            # Extraer la ventana de procesos
            window = flux_values[i : i + scale]
            
            # Sumar los vectores de proceso en la ventana para obtener el "proceso combinado"
            v_combined = np.sum(window, axis=0)
            
            # Clasificar el proceso combinado
            classifications = classify_process(v_combined, S)
            
            # Contar si pertenece al Dominio Cognitivo
            if "Cognitive Domain" in classifications:
                cognitive_domain_count += 1
        
        results[scale] = cognitive_domain_count
        print(f"Escala {scale}: {cognitive_domain_count} procesos de Dominio Cognitivo encontrados.")

    # --- Visualización de los Resultados ---
    scales = list(results.keys())
    counts = list(results.values())

    # Encontrar la escala con el máximo conteo
    optimal_scale = max(results, key=results.get)
    max_count = results[optimal_scale]

    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(12, 7))
    plt.plot(scales, counts, marker='o', linestyle='-', color='b', label='Conteo de Dominio Cognitivo')
    
    # Resaltar el punto óptimo
    plt.axvline(x=optimal_scale, color='r', linestyle='--', label=f'Escala Óptima = {optimal_scale} (Conteo: {max_count})')
    plt.scatter(optimal_scale, max_count, color='red', s=100, zorder=5)

    plt.title('Análisis de la Escala de Tiempo de Estabilidad', fontsize=16, fontweight='bold')
    plt.xlabel('Escala de Tiempo (Tamaño de la Ventana)', fontsize=12)
    plt.ylabel('Frecuencia de Procesos de "Dominio Cognitivo"', fontsize=12)
    plt.xticks(scales)
    plt.legend()

    # --- Guardar la figura ---
    out_dir="./visualizations/process_classification" 
    filename="time_scale_analysis.png"
    os.makedirs(out_dir, exist_ok=True)
    save_path = os.path.join(out_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Figura guardada en: {save_path}")

    plt.show()

    return results

################################################################################################# 
# Función para crear series de tiempo agrupadas  
def new_time_series(time_series, flux_vector, scale):
    """
    Agrupa los resultados de una simulación en una escala de tiempo mayor.

    Args:
        time_series (pd.DataFrame): DataFrame de concentraciones de la simulación.
        flux_vector (pd.DataFrame): DataFrame de flujos de la simulación.
        scale (int): La escala o tamaño de la ventana para agrupar.

    Returns:
        tuple: Una tupla conteniendo (ts_coarse_df, fv_coarse_df), los nuevos
               DataFrames agrupados.
    """
    if scale < 1:
        return time_series, flux_vector

    # Extraer valores como arrays de NumPy para eficiencia
    ts_values = time_series.to_numpy()
    fv_values = flux_vector.iloc[:, 1:].to_numpy() # Excluir columna de tiempo

    new_ts_rows = []
    new_fv_rows = []

    # Iterar en bloques no superpuestos del tamaño de la escala
    for i in range(0, len(ts_values) - scale + 1, scale):
        # La ventana de flujos a sumar
        flux_window = fv_values[i : i + scale]
        
        # El proceso combinado es la suma de los flujos en la ventana
        combined_flux = np.sum(flux_window, axis=0)
        
        # El punto de tiempo y las concentraciones corresponden al final de la ventana
        end_of_window_ts = ts_values[i + scale - 1]
        
        # Añadir la nueva fila de flujo (tiempo + flujos sumados)
        new_fv_rows.append(np.insert(combined_flux, 0, end_of_window_ts[0]))
        
        # Añadir la nueva fila de concentraciones
        new_ts_rows.append(end_of_window_ts)

    # Crear nuevos DataFrames con las columnas originales
    ts_coarse_df = pd.DataFrame(new_ts_rows, columns=time_series.columns)
    fv_coarse_df = pd.DataFrame(new_fv_rows, columns=flux_vector.columns)
    
    return ts_coarse_df, fv_coarse_df

def find_min_cognitive_interval(rn, x0, spec_vector, rate_list='mak', t_max_values=None, n_steps_values=None):
    from pyCOT.simulations import simulation
    from pyCOT.plot_dynamics import plot_flux_with_domain_intervals
    
    if t_max_values is None:
        t_max_values = [50, 75, 100, 125]  # puedes modificar
    if n_steps_values is None:
        n_steps_values = [500, 1000, 2000, 3000]  # puedes modificar
    
    best_interval = float('inf')
    best_config = None
    best_flux_vector = None
    best_time_series = None
    
    for t_max in t_max_values:
        for n_steps in n_steps_values:
            # Ejecutar simulación
            time_series, flux_vector = simulation(
                rn, rate=rate_list, spec_vector=spec_vector, x0=x0,
                t_span=(0, t_max), n_steps=n_steps+1
            )
            
            # Obtener intervalos del Cognitive Domain
            intervals = plot_flux_with_domain_intervals(
                flux_vector, rn.stoichiometry_matrix(),
                title=f"t_max={t_max}, n_steps={n_steps}", 
                save_figure=False, show_fig=False
            )
            
            # intervals devuelve lista de tuplas: [(start, end), ...]
            if intervals:
                min_interval = min([end-start for start, end in intervals])
                if min_interval < best_interval:
                    best_interval = min_interval
                    best_config = (t_max, n_steps)
                    best_flux_vector = flux_vector
                    best_time_series = time_series
    
    print(f"Mejor configuración: t_max={best_config[0]}, n_steps={best_config[1]}")
    print(f"Intervalo mínimo encontrado: {best_interval}")
    
    return best_time_series, best_flux_vector, best_config, best_interval

import numpy as np
from pyCOT.simulations import simulation

def find_natural_time_scale(rn, x0=None, spec_vector=None, rate_list='mak', 
                            t_max_values=None, n_steps_values=None, tol=1e-4):
    """
    Busca la escala de tiempo natural (t_max, n_steps) de la simulación que genere el 
    mínimo periodo T.

    Args:
        rn: Red de reacciones (ReactionNetwork object)
        x0: Condiciones iniciales
        spec_vector: Vector de especificidad de reacciones
        rate_list: Cinéticas
        t_max_values: Lista de tiempos máximos a probar
        n_steps_values: Lista de n_steps a probar
        tol: Tolerancia para la clasificación de procesos

    Returns:
        best_time_series, best_flux_vector, best_config, min_interval
    """
    
    S = rn.stoichiometry_matrix()
    
    if t_max_values is None:
        t_max_values = [50, 75, 100, 125]
    if n_steps_values is None:
        n_steps_values = [500, 1000, 2000, 3000]

    best_interval = float('inf') # Se inicializa con infinito para poder comparar y guardar el menor intervalo
    best_config = None
    best_flux_vector = None
    best_time_series = None

    # Iterar sobre todas las combinaciones de t_max y n_steps
    for t_max in t_max_values:
        for n_steps in n_steps_values:
            # Ejecutar simulación con la configuración actual
            time_series, flux_vector = simulation(
                rn, rate=rate_list, spec_vector=spec_vector, x0=x0,
                t_span=(0, t_max), n_steps=n_steps+1
            )

            # Clasificar cada vector de flujo para cada instante de tiempo
            cognitive_times = []
            for i in range(len(flux_vector)):
                v = flux_vector.iloc[i, 1:].values  # excluyendo columna Time
                classes = classify_process(v, S, tol=tol)
                if "Cognitive Domain" in classes:
                    cognitive_times.append(flux_vector.iloc[i, 0])  # Guardar el tiempo

            # Detectar intervalos consecutivos del Cognitive Domain
            if cognitive_times:
                cognitive_times = np.array(cognitive_times) # Convertir a array de NumPy
                diffs = np.diff(cognitive_times) # Diferencias entre tiempos consecutivos

                # Considerar un salto mayor a dt como fin de un intervalo
                # dt = cognitive_times[1] - cognitive_times[0] # Paso de tiempo típico entre dos instantes consecutivos que fueron clasificados dentro del Cognitive Domain
                dt = cognitive_times[1] - cognitive_times[0] if len(cognitive_times) > 1 else 0
                splits = np.where(diffs > 1.5*dt)[0]  # 1.5*dt para tolerancia. Si la diferencia es un poco más grande que dt (por ruido numérico), no corta el intervalo.
                
                intervals = []
                start_idx = 0
                for split in splits:
                    intervals.append((cognitive_times[start_idx], cognitive_times[split]))
                    start_idx = split+1
                # Se guardan todos los intervalos como (t_start, t_end).    
                intervals.append((cognitive_times[start_idx], cognitive_times[-1]))  # último intervalo

                # Intervalo mínimo en esta configuración
                min_interval_current = min([end-start for start, end in intervals])
                if min_interval_current < best_interval:
                    best_interval = min_interval_current
                    best_config = (t_max, n_steps)
                    best_flux_vector = flux_vector
                    best_time_series = time_series

    if best_config is not None:
        print(f"Mejor configuración: t_max={best_config[0]}, n_steps={best_config[1]}")
        return best_time_series, best_flux_vector, best_config, best_interval
    else:
        print("⚠️ No se encontró ninguna configuración válida para la escala de tiempo natural.")
        return None, None, None, None

    # print(f"Mejor configuración: t_max={best_config[0]}, n_steps={best_config[1]}")
    # print(f"Intervalo mínimo del Cognitive Domain: {best_interval:.5f}")

    # return best_time_series, best_flux_vector, best_config, best_interval
