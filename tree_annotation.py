import json
import random

import numpy as np
import pandas as pd
import anndata as ad

import math
import tkinter as tk
from bokeh.models import canvas


def tree_to_umap(umap_file):
    cell_type_to_umap = {}
    adata = ad.read_h5ad(umap_file)
    level_columns = [col for col in adata.obs.columns if col.startswith('level_')]

    for level in level_columns:
        unique_cell_types = adata.obs[level].unique()
        for cell_type in unique_cell_types:
            bool_indices = adata.obs[level] == cell_type
            int_indices = np.where(bool_indices)[0]
            umap_values = adata.obsm['X_umap'][int_indices]

            if cell_type not in cell_type_to_umap:
                cell_type_to_umap[cell_type] = []

            cell_type_to_umap[cell_type].extend(umap_values.tolist())
    return cell_type_to_umap


def get_node_levels(umap_file):
    global node
    node_levels = {}
    adata = ad.read_h5ad(umap_file)
    level_columns = [col for col in adata.obs.columns if col.startswith('level_')]

    for level in level_columns:
        unique_cell_types = adata.obs[level].unique()

        key = ','.join(map(str, unique_cell_types))
        node_levels[key] = level.split('_')[1]
    for key, value in node_levels.copy().items():
        if len(key) > 1:
            for cell_type in key.split(','):
                node_levels[cell_type] = value
            del node_levels[key]
    return node_levels


def sort_node(selected_nodes):
    global node_levels
    try:
        sorted_nodes = sorted(
            selected_nodes,
            key=lambda node: str(node_levels.get(node, ''))
        )
    except ValueError as e:
        print(f"Error converting levels to float: {e}")
        return selected_nodes
    print(sorted_nodes)
    return sorted_nodes



def get_node_positions(position_file, mapping_file):
    with open(position_file) as f:
        node_positions = json.load(f)
    with open(mapping_file) as w:
        cell_type_map = json.load(w)
    node_position = {}
    for key, value in node_positions.items():
        cell_type = cell_type_map.get(key, '')
        if cell_type:
            node_position[cell_type] = value
        else:
            raise KeyError(f"Key '{key}' in node_positions not found in cell_type_map.")
    return node_position

def create_circle_button(parent, x, y, radius, name, command=None):
    if name not in node_colors.keys():
        node_color = generate_random_color()
        node_colors[name] = node_color
    else:
        node_color = node_colors[name]
    button = tk.Button(
        parent,
        text='',
        command=lambda: command(name),
        font=("Arial", 5),
        bg=node_color,
        fg="white",
        relief="flat",
        bd=0
    )
    button.place(
        x=0.85 * x - radius + 40,
        y=0.85 * y - radius + 40,
        width=2 * radius,
        height=2 * radius
    )

    parent.update_idletasks()
    button_x = button.winfo_x()
    button_y = button.winfo_y()
    button_width = button.winfo_width()

    label = tk.Label(
        parent,
        text=name,
        font=("Arial", 10),
        bg="grey",
        fg="white"
    )
    label.place(
        x=button_x + button_width / 2 - label.winfo_reqwidth() / 2,
        y=button_y - 20
    )
    return button

def draw_arrow(canvas, x1, y1, x2, y2, color="white"):

    canvas.create_line(x1, y1, x2, y2, arrow=tk.LAST, fill=color, width=2)

def create_tree_button():
    global canvas
    for node_name, (x, y) in node_positions.items():
        tkinter_x, tkinter_y = convert_coordinates(x, y, min_x, max_x, min_y, max_y, gui_width, gui_height)
        create_circle_button(
            root, tkinter_x, tkinter_y, radius=6, name=node_name,
            command=lambda name=node_name: button_on_click(None, name, text_widget)
        )

def creat_background():
    global canvas
    for _, row in edge.iterrows():
        parent = row['parents']
        subnode = row['subnodes']
        if parent in node_positions and subnode in node_positions:
            parent_pos = node_positions[parent]
            subnode_pos = node_positions[subnode]

            x1, y1 = convert_coordinates(parent_pos[0], parent_pos[1], min_x, max_x, min_y, max_y, gui_width,
                                         gui_height)
            x2, y2 = convert_coordinates(subnode_pos[0], subnode_pos[1], min_x, max_x, min_y, max_y, gui_width,
                                         gui_height)
            x1 = 0.85 * x1 + 36,
            y1 = 0.85 * y1 + 42,
            x2 = 0.85 * x2 + 36,
            y2 = 0.85 * y2 + 42,
            draw_arrow(canvas, x1, y1, x2, y2)


def create_continue_button():
    bottom_frame = tk.Frame(root, bg="grey", height=40)
    bottom_frame.place(relx=0, rely=1.0, anchor="sw", width=gui_width+5)

    continue_button = tk.Button(bottom_frame, text="Continue", command=on_continue, width=15, bg="grey", fg="white")
    continue_button.pack(pady=5)


def create_back_button():
    bottom_frame_back = tk.Frame(root, bg="grey")
    bottom_frame_back.place(relx=0, rely=1.0, anchor="sw", width=gui_width+5)

    back_button = tk.Button(bottom_frame_back, text="Back", command=on_back, width=15, bg="white", fg="grey")
    back_button.pack(pady=10)

def convert_coordinates(x, y, min_x, max_x, min_y, max_y, gui_width, gui_height):
    tkinter_x = (x - min_x) / (max_x - min_x) * gui_width
    tkinter_y = gui_height - ((y - min_y) / (max_y - min_y) * gui_height)
    return tkinter_x, tkinter_y


def create_info(parent):

    frame = tk.Frame(parent, bg="lightgray", width=300, height=150)
    frame.place(relx=1.0, rely=1.0, anchor="se")

    label = tk.Label(frame, text="Selected Cell Type:", bg="lightgray", font=("Arial", 10, "bold"))
    label.pack(anchor="nw", padx=5, pady=5)

    text_widget = tk.Text(frame, height=8, width=20, bg="#f9f9f9", wrap="none", font=("Arial", 9))
    text_widget.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    scrollbar = tk.Scrollbar(frame, command=text_widget.yview)
    text_widget.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    return text_widget

def create_result(parent):
    frame = tk.Frame(parent, bg="white", width=300, height=150)
    frame.place(relx=1.0, rely=1.0, anchor="se")

    result = tk.Text(frame, height=8, width=20, bg="white", wrap="none", font=("Arial", 9))
    result.pack(fill=tk.BOTH, expand=True, padx=5, pady=18)

    scrollbar = tk.Scrollbar(frame, command=result.yview)
    result.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    return result

selected_nodes = []
node_colors = {}

global_min_x, global_max_x, global_min_y, global_max_y = None, None, None, None

def calculate_global_bounds(selected_name):

    all_coords = [coord for node_name in selected_name for coord in cell_type_to_umap[node_name]]
    min_x = min(all_coords, key=lambda x: x[0])[0]
    max_x = max(all_coords, key=lambda x: x[0])[0]
    min_y = min(all_coords, key=lambda x: x[1])[1]
    max_y = max(all_coords, key=lambda x: x[1])[1]
    return min_x, max_x, min_y, max_y

def initialize_scaling(selected_name, scale_factor=1.5):

    global global_min_x, global_max_x, global_min_y, global_max_y
    global canvas_width, canvas_height
    global node

    if selected_name:
        min_x, max_x, min_y, max_y = calculate_global_bounds(node)

        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        range_x = (max_x - min_x) * scale_factor
        range_y = (max_y - min_y) * scale_factor
        global_min_x = center_x - range_x / 2
        global_max_x = center_x + range_x / 2
        global_min_y = center_y - range_y / 2
        global_max_y = center_y + range_y / 2

def convert_umap(node_name):
    global global_min_x, global_max_x, global_min_y, global_max_y

    if None in (global_min_x, global_max_x, global_min_y, global_max_y):
        raise ValueError("Global bounds have not been initialized. Call `initialize_scaling` first.")

    scaled_x = [
        (x - global_min_x) / (global_max_x - global_min_x) * canvas_width
        for x, y in cell_type_to_umap[node_name]
    ]
    scaled_y = [
        (y - global_min_y) / (global_max_y - global_min_y) * canvas_height
        for x, y in cell_type_to_umap[node_name]
    ]
    scaled_coord = [(new_x, new_y) for new_x, new_y in zip(scaled_x, scaled_y)]
    return scaled_coord


def generate_random_color():
    return f"#{random.randint(0, 255):02x}{random.randint(0, 255):02x}{random.randint(0, 255):02x}"

def button_on_click(event, node_name, text_widget):
    global selected_nodes
    global node_colors
    if node_name not in selected_nodes:
        selected_nodes.append(node_name)
        text_widget.insert(tk.END, f"{node_name}\n")

    else:
        selected_nodes.remove(node_name)

        lines = text_widget.get("1.0", tk.END).splitlines()
        for i in range(len(lines) - 1, -1, -1):
            if lines[i] == node_name:
                text_widget.delete(f"{i+1}.0", f"{i+1}.0 lineend+1c")
                break


def highlight_umap(canvas, scaled_coord, node_color):
    for coord in scaled_coord:
        if len(coord) != 2:
            print(f"Invalid coordinate format: {coord}")
            continue

        x, y = coord
        if isinstance(x, (int, float)) and isinstance(y, (int, float)):
            canvas.create_oval(
                x - 5, y - 5, x + 5, y + 5,
                fill=node_color,
                outline=node_color
            )
        else:
            print(f"Invalid coordinate value: {coord}")

def on_continue():
    global selected_nodes, node_colors, canvas_width, canvas_height
    global node, node_levels
    for widget in root.winfo_children():
        widget.destroy()

    root.configure(bg='white')

    bottom_frame_back = tk.Frame(root, bg="white")
    bottom_frame_back.pack(side="bottom", fill="x", pady=10)
    back_button = tk.Button(bottom_frame_back, text="Back", command=on_back, width=15, bg="white", fg="grey")
    back_button.pack(pady=10)

    umap_canvas = tk.Canvas(root, width=gui_width, height=gui_height, bg="white")
    umap_canvas.pack(fill="both", expand=True)
    umap_canvas.update_idletasks()
    canvas_width = umap_canvas.winfo_width()
    canvas_height = umap_canvas.winfo_height()
    initialize_scaling(node)
    result = create_result(root)
    if not selected_nodes:

        for node_name in node:
            scaled_coord = convert_umap(node_name)
            highlight_umap(umap_canvas, scaled_coord, 'lightgrey')
    else:
        remain_node = list(set(node) - set(selected_nodes))
        sorted_node = sort_node(selected_nodes)
        for node_name in remain_node:
            scaled_coord = convert_umap(node_name)
            highlight_umap(umap_canvas, scaled_coord, 'lightgrey')

        for node_name in sorted_node:
            node_color = node_colors.get(node_name, "black")
            print(f"Node: {node_name}, Color: {node_color}")
            scaled_coord = convert_umap(node_name)
            highlight_umap(umap_canvas, scaled_coord, node_color)
        print("...........")

        for idx, node_name in enumerate(sorted_node):
            node_name = node_name.strip()
            node_color = node_colors.get(node_name, "black").strip()

            print(f"Node: {node_name}, Color: {node_color}")

            tag_name = f"{node_name}_tag_{idx}"
            result.tag_config(tag_name, foreground=node_color)
            print(f"{tag_name}: {node_color}")

            result.insert(tk.END, f"{node_name} ")

            start_index = result.index(tk.END)
            result.insert(tk.END, "●")
            end_index = result.index(tk.END)
            result.tag_add(tag_name, start_index + "-1c", end_index)

            result.insert(tk.END, '\n')


def on_back():
    global selected_nodes, text_widget, canvas

    for widget in root.winfo_children():
        widget.destroy()
    selected_nodes = []
    root.configure(bg='grey')
    canvas = tk.Canvas(root, width=gui_width, height=gui_height, bg="grey")
    canvas.pack(fill="both", expand=True)
    create_tree_button()
    creat_background()
    create_continue_button()
    text_widget = create_info(root)
    return text_widget


umap_file = r'dataset\catree_data.h5ad'
position_file = r'dataset\node_positions.json'
mapping_file = r'dataset\cell_type.json'
edge_file = r'dataset\edge.csv'

cell_type_to_umap = tree_to_umap(umap_file=umap_file)
node_position = get_node_positions(position_file=position_file, mapping_file=mapping_file)
node_levels = get_node_levels(umap_file)

edge = pd.read_csv(edge_file)
node = list(node_position.keys())
coordinates = cell_type_to_umap["cell"]
x = [coord[0] for coord in coordinates]
y = [coord[1] for coord in coordinates]
min_x, max_x = min(x), max(x)
min_y, max_y = min(y), max(y)

node_positions = {}
theta = math.pi / 2
rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                            [np.sin(theta), np.cos(theta)]])
for node_name, (x, y) in node_position.items():
    new_position = rotation_matrix @ np.array([x, y])
    node_positions[node_name] = (new_position[0], new_position[1])


x_coords = [pos[0] for pos in node_positions.values()]
y_coords = [pos[1] for pos in node_positions.values()]
min_x, max_x = min(x_coords), max(x_coords)
min_y, max_y = min(y_coords), max(y_coords)


gui_width = 1400
gui_height = 700

root = tk.Tk()
root.title("GUI")
root.geometry(f"{gui_width}x{gui_height}")
root.configure(bg="grey")

canvas = tk.Canvas(root, width=gui_width, height=gui_height, bg="grey")
canvas.pack(fill="both", expand=True)
creat_background()
create_tree_button()
create_continue_button()
text_widget = create_info(root)
root.mainloop()

