# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 14:14:30 2025

@author: jannatul
FLARE version 1.1
"""


import pandas as pd
import rasterio
from rasterio.transform import rowcol
from shapely.geometry import Point
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import numpy as np
import random
import tkinter as tk
import pyproj
from pyproj import transformer

# Initialize tkinter
%matplotlib tk
root = tk.Tk()
root.withdraw()

# Input file paths
LANDCOVER_FILE = "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/Friesland_2025.tif"
BIRD_MOVEMENT_FILE = "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/Filtered_Vanellus_V.csv"

# BirdAgent class
class BirdAgent:
    def __init__(self, unique_id, movement_data):
        self.unique_id = unique_id
        self.movement_data = movement_data  # Bird's recorded movement data
        self.current_step = 0  # Track movement data index
        self.current_position = None
        self.previous_position = None
        self.habitat_preference = {}
        self.behavior_log = []  # To record each step's habitat and behavior

    def step(self, model, sensitivity=None):
        """Move the bird based on movement data and record habitat preferences."""
        if self.current_step < len(self.movement_data):
            # Read the bird's location
            data = self.movement_data.iloc[self.current_step]
            timestamp = data['timestamp']
            current_position = Point(data['location.long'], data['location.lat'])
    
            # Get the habitat type at the current position
            habitat_value = model.get_habitat_at_position(current_position)
    
            # Adjust step length based on sensitivity
            step_length_multiplier = 1.0  # Default multiplier used if there is no sensitivity score for a habitat, 1=neutral
            if habitat_value is not None and sensitivity is not None:
                habitat_sensitivity = sensitivity.get(habitat_value, 1.0)
                step_length_multiplier = habitat_sensitivity 
                print(f"Bird {self.unique_id}: Habitat {habitat_value}, Sensitivity {habitat_sensitivity}, Step Multiplier {step_length_multiplier}")
    
            # Calculate step length and turning angle
            step_length = 0
            turning_angle = None

            if self.previous_position:
                dx = (current_position.x - self.previous_position.x) * step_length_multiplier
                dy = (current_position.y - self.previous_position.y) * step_length_multiplier
                step_length = np.sqrt(dx**2 + dy**2)

                if self.current_position:
                    prev_dx = self.current_position.x - self.previous_position.x
                    prev_dy = self.current_position.y - self.previous_position.y
            
                    dot_product = prev_dx * dx + prev_dy * dy
                    magnitude_a = np.sqrt(prev_dx**2 + prev_dy**2)
                    magnitude_b = np.sqrt(dx**2 + dy**2)

                    if magnitude_a > 0 and magnitude_b > 0:
                        cosine_angle = dot_product / (magnitude_a * magnitude_b)
                        cosine_angle = np.clip(cosine_angle, -1.0, 1.0)  # <--- safe clipping
                        turning_angle = np.arccos(cosine_angle)
                        turning_angle = np.degrees(turning_angle)

    
            # Update positions
            self.previous_position = self.current_position
            self.current_position = current_position
    
            # Record habitat preference and behavior
            if habitat_value is not None:
                self.habitat_preference[habitat_value] = self.habitat_preference.get(habitat_value, 0) + 1
    
            self.behavior_log.append({
                "Step": self.current_step + 1,
                "Timestamp": timestamp,
                "Habitat Class": habitat_value,
                "Longitude": current_position.x,
                "Latitude": current_position.y,
                "Step Length": step_length,
                "Turning Angle": turning_angle
            })
    
            self.current_step += 1

#SmartBirdAgent class for hypothetical model
class SmartBirdAgent(BirdAgent):
    """BirdAgent that moves based on habitat sensitivity."""
    
    def step(self, model, sensitivity=None):
        """Move the bird intelligently based on habitat preferences."""
        if self.current_step == 0:
            # Initialize at the real first position
            data = self.movement_data.iloc[self.current_step]
            current_position = Point(data['location.long'], data['location.lat'])
            self.current_position = current_position
            self.previous_position = None
            self.current_step += 1
            return

        if self.current_position is None:
            return
        
        directions = [
            (0.0005, 0), (0.0005, 0.0005), (0, 0.0005), (-0.0005, 0.0005),
            (-0.0005, 0), (-0.0005, -0.0005), (0, -0.0005), (0.0005, -0.0005)
        ]##The direction vector string indicating 8 directions
        
        candidates = []
        for dx, dy in directions:
            candidate_pos = Point(self.current_position.x + dx, self.current_position.y + dy)
            habitat_value = model.get_habitat_at_position(candidate_pos)
        
            if habitat_value is not None:
                # if normal valid move record location
                score = sensitivity.get(habitat_value, 1.0)
                candidates.append((score, candidate_pos, habitat_value))
            else:
                # Outside raster: Bounce back within the map boundary
                bounce_pos = Point(self.current_position.x - dx, self.current_position.y - dy)
                bounce_habitat_value = model.get_habitat_at_position(bounce_pos)
        
                if bounce_habitat_value is not None:
                    bounce_score = sensitivity.get(bounce_habitat_value, 1.0)
                    candidates.append((bounce_score, bounce_pos, bounce_habitat_value))
        
        # Now choose among candidates
        if not candidates:
            # No valid moves even after bouncing
            self.behavior_log.append({
                "Step": self.current_step + 1,
                "Timestamp": None,
                "Habitat Class": None,
                "Longitude": self.current_position.x,
                "Latitude": self.current_position.y,
                "Step Length": 0,
                "Turning Angle": None
            })
            self.current_step += 1
            return


        # Correct sorting - sort by score descending
        candidates.sort(key=lambda x: x[0], reverse=True)

        # pick randomly among top 26 
        top_n = min(6, len(candidates))
        best_score, best_position, best_habitat = random.choice(candidates[:top_n])
        
        
        self.previous_position = self.current_position
        self.current_position = best_position
        
        step_length = np.sqrt((self.current_position.x - self.previous_position.x)**2 +
                              (self.current_position.y - self.previous_position.y)**2)
        
        turning_angle = None
        
        self.habitat_preference[best_habitat] = self.habitat_preference.get(best_habitat, 0) + 1
        
        self.behavior_log.append({
            "Step": self.current_step + 1,
            "Timestamp": None, #self.movement_data.iloc[self.current_step]['timestamp'],
            "Habitat Class": best_habitat,
            "Longitude": self.current_position.x,
            "Latitude": self.current_position.y,
            "Step Length": step_length,
            "Turning Angle": turning_angle
        })
        
        self.current_step += 1


# BirdHabitatModel class
class BirdHabitatModel:
    def __init__(self, landcover_file, bird_movement_file):
        self.landcover = rasterio.open(landcover_file)
        self.landcover_data = self.landcover.read(1)
        # Merge Ocean (0) into Salt Water (17) (Case specific decision)
        ##To be used when there's more than one class that are similar
        #and merging them makes sence based on the species behaviour
        #Mute if no merge require
        self.landcover_data[self.landcover_data == 0] = 17

        # Initialize agents
        self.bird_data = pd.read_csv(bird_movement_file)
        self.bird_agents = {}
        self.bird_colors = {}

        # Assign colors to birds for the python animation plot
        color_list = ['red', 'pink', 'orange', 'purple']
        for idx, bird_id in enumerate(self.bird_data['tag.local.identifier'].unique()):
            movement_data = self.bird_data[self.bird_data['tag.local.identifier'] == bird_id]
            agent = BirdAgent(bird_id, movement_data)
            self.bird_agents[bird_id] = agent
            self.bird_colors[bird_id] = color_list[idx % len(color_list)]

    def get_habitat_at_position(self, position):
        """Get the habitat class from the raster at the given position."""
        lon, lat = position.x, position.y
        try:
            row, col = rowcol(self.landcover.transform, lon, lat)
            habitat_value = self.landcover_data[row, col]
            return habitat_value
        except IndexError:
            print(f"Position {position} is outside the raster bounds.")
            return None

    def step(self):
        """Advance the model by one time step."""
        for agent in self.bird_agents.values():
            agent.step(self, sensitivity=self.sensitivity if hasattr(self, 'sensitivity') else None)


    def get_results(self):
        """Retrieve habitat preference results and behavior logs."""
        results = []
        for bird_id, agent in self.bird_agents.items():
            for log in agent.behavior_log:
                results.append({
                    "Bird ID": bird_id,
                    **log  # Include all log fields
                })
        return results

   
    def calculate_sensitivity(self):
        """Calculate sensitivity for each landcover type."""
        # Count usage for each habitat type
        habitat_usage = {}
        for agent in self.bird_agents.values():
            for habitat, count in agent.habitat_preference.items():
                habitat_usage[habitat] = habitat_usage.get(habitat, 0) + count
    
        # Count availability of each habitat type in the landcover map
        unique, counts = np.unique(self.landcover_data, return_counts=True)
        habitat_availability = dict(zip(unique, counts))
    
        # Calculate sensitivity scores
        sensitivity = {}
        total_usage = sum(list(habitat_usage.values()))  
        total_availability = self.landcover_data.size
        for habitat, usage_count in habitat_usage.items():
            usage_proportion = usage_count / total_usage if total_usage > 0 else 0
            availability_proportion = habitat_availability.get(habitat, 0) / total_availability if total_availability > 0 else 0
            sensitivity[habitat] = usage_proportion / availability_proportion if availability_proportion > 0 else 0
    
        print("Sensitivity Scores:", sensitivity)
        return sensitivity



# HypotheticalBirdHabitatModel
class HypotheticalBirdHabitatModel(BirdHabitatModel):
    def __init__(self, landcover_file, bird_movement_file, sensitivity):
        super().__init__(landcover_file, bird_movement_file)
        self.sensitivity = sensitivity  

        # Modify the landcover map
        most_used_habitat = max(sensitivity, key=sensitivity.get)
        self.hypothetical_landcover_data = self.landcover_data.copy()
        
        class_indices = np.argwhere(self.hypothetical_landcover_data == most_used_habitat)
        num_to_convert = int(len(class_indices) * 0.33)  # Replace % for landuse modification
        
        selected_indices = random.sample(list(class_indices), num_to_convert)
        for row, col in selected_indices:
            self.hypothetical_landcover_data[row, col] = 12  # New habitat class (to which the preferred one will be converted)
        
        print(f"Replaced 39% of class {most_used_habitat} with class 26.")

        # Replace bird agents with SmartBirdAgents
        self.bird_agents = {}
        self.bird_colors = {}

        color_list = ['red', 'pink', 'orange', 'purple']
        
        for idx, bird_id in enumerate(self.bird_data['tag.local.identifier'].unique()):
            movement_data = self.bird_data[self.bird_data['tag.local.identifier'] == bird_id]
            agent = SmartBirdAgent(bird_id, movement_data)  # <--- Use SmartBirdAgent here!
            self.bird_agents[bird_id] = agent
            self.bird_colors[bird_id] = color_list[idx % len(color_list)]

    def get_habitat_at_position(self, position):
        """Get the habitat class from the hypothetical raster at the given position."""
        lon, lat = position.x, position.y
        try:
            row, col = rowcol(self.landcover.transform, lon, lat)
            habitat_value = self.hypothetical_landcover_data[row, col]
            return habitat_value
        except IndexError:
            print(f"Position {position} is outside the raster bounds.")
            return None


# Combined Animation Function: showing both scenario side by side in one animation panel
def animate(i):
    """Update the animation for frame i."""
    current_model.step()
    hypothetical_model.step()

    
    ax1.clear()
    ax2.clear()

    # Plot current scenario
    ax1.imshow(current_model.landcover_data, cmap='nipy_spectral', norm=norm_current)
    ax1.set_title(f"Current Scenario: Time Step {i + 1}")
    for bird_id, agent in current_model.bird_agents.items():
        if agent.current_position:
            lon, lat = agent.current_position.x, agent.current_position.y
            row, col = rowcol(current_model.landcover.transform, lon, lat)
            ax1.scatter(col, row, color=current_model.bird_colors[bird_id], label=f"Bird {bird_id}", s=80)
    ax1.legend()

    # Plot hypothetical scenario
    ax2.imshow(hypothetical_model.hypothetical_landcover_data, cmap='nipy_spectral', norm=norm_hypothetical)
    ax2.set_title(f"Hypothetical Scenario: Time Step {i + 1}")
    for bird_id, agent in hypothetical_model.bird_agents.items():
        if agent.current_position:
            lon, lat = agent.current_position.x, agent.current_position.y
            row, col = rowcol(hypothetical_model.landcover.transform, lon, lat)
            ax2.scatter(col, row, color=hypothetical_model.bird_colors[bird_id], label=f"Bird {bird_id}", s=80)
    ax2.legend()


# Main Function
if __name__ == "__main__":
    # Initialize the current model
    current_model = BirdHabitatModel(LANDCOVER_FILE, BIRD_MOVEMENT_FILE)

    # Run the current model
    num_time_steps = 30000 ##To be adjusted considering available data and the length of period of life cycle that covers
    for step in range(num_time_steps):
        current_model.step()

    # Calculate sensitivity
    sensitivity = current_model.calculate_sensitivity()

    # Export sensitivity scores
    sensitivity_df = pd.DataFrame(list(sensitivity.items()), columns=["Habitat Class", "Sensitivity Score"])
    sensitivity_df = sensitivity_df.sort_values(by="Sensitivity Score", ascending=False)
    sensitivity_df.to_excel("D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Sensitivity_current_41to12.xlsx", index=False)
    print("Sensitivity scores exported successfully!")

    # Initialize the hypothetical model with sensitivity
    hypothetical_model = HypotheticalBirdHabitatModel(LANDCOVER_FILE, BIRD_MOVEMENT_FILE, sensitivity)

    # Run the hypothetical model
    for step in range(num_time_steps):
        hypothetical_model.step()
    # Calculate new sensitivity for hypothetical scenario
    hypothetical_sensitivity = hypothetical_model.calculate_sensitivity()

    # Export hypothetical sensitivity
    hypothetical_sensitivity_df = pd.DataFrame(list(hypothetical_sensitivity.items()), columns=["Habitat Class", "Sensitivity Score"])
    hypothetical_sensitivity_df = hypothetical_sensitivity_df.sort_values(by="Sensitivity Score", ascending=False)
    hypothetical_sensitivity_df.to_excel("D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Sensitivity_Hypo_41to12.xlsx", index=False)
    print("Hypothetical sensitivity scores exported successfully!")

    # Export results
    current_df = pd.DataFrame(current_model.get_results())
    hypothetical_df = pd.DataFrame(hypothetical_model.get_results())
    
    current_df.to_excel("D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Current_41to12.xlsx", index=False)
    hypothetical_df.to_excel("D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Hypo_41to12.xlsx", index=False)
    print("Results exported successfully!")
####Turn the comment off in case you want animation view in Python 
###Muted as it requires better computing capacity
""""
    # Set up subplots for animation
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    norm_current = Normalize(vmin=current_model.landcover_data.min(), vmax=current_model.landcover_data.max())
    norm_hypothetical = Normalize(vmin=hypothetical_model.hypothetical_landcover_data.min(),
                                  vmax=hypothetical_model.hypothetical_landcover_data.max())
    
    # Run animation
    ani = animation.FuncAnimation(fig, animate, frames=100, repeat=False)
    # Save as MP4
    output_path = "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/combined_animation.mp4"
    ani.save(output_path, writer='ffmpeg', fps=5)

    print("Animation saved as MP4!")
    plt.tight_layout()
    plt.show()
"""""

