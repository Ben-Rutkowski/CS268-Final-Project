# importing libraries
import pandas as pd
from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry

cities = pd.read_csv(r'CS268_EVDemand.csv')

# converting Longitude and Latitude to xy coordinates
geometry = [Point(xy) for xy in zip(cities['Longitude'], cities['Latitude'])]

# Translate coordinates to positive values
x = [geometry[i].x for i in range(len(geometry))]
x = x-min(x)
y = [geometry[i].y for i in range(len(geometry))]
y = y-min(y)

cities.insert(4, 'x', x)
cities.insert(5, 'y', y)

# save new csv file
cities.to_csv(r'CS268_EVDemand_xy_normalized.csv', encoding='utf-8', header='true')