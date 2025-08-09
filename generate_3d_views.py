#!/usr/bin/env python3
"""
3D Model View Generator for Rhino .3dm files
Generates orthographic and isometric wireframe views from actual 3DM geometry
"""

import matplotlib
matplotlib.use('Agg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse
from pathlib import Path

def find_3dm_file(directory='.'):
    """Find the first .3dm file in the specified directory"""
    for file in Path(directory).glob('*.3dm'):
        return str(file)
    return None

def sample_curve_points(curve, count=50):
    """Sample points along a curve"""
    points = []
    try:
        # Try multiple approaches for different curve types
        if hasattr(curve, 'Domain'):
            domain = curve.Domain
            # Fix: Use T0 and T1 instead of Min and Max for rhino3dm intervals
            t_start = domain.T0
            t_end = domain.T1
            for i in range(count):
                t = t_start + (t_end - t_start) * i / (count - 1)
                if hasattr(curve, 'PointAt'):
                    pt = curve.PointAt(t)
                    if pt:
                        points.append([pt.X, pt.Y, pt.Z])
                        
        elif hasattr(curve, 'Points') and hasattr(curve.Points, 'Count'):
            # NURBS curve with control points
            point_count = curve.Points.Count
            step = max(1, point_count // count)
            for i in range(0, point_count, step):
                pt = curve.Points[i].Location
                points.append([pt.X, pt.Y, pt.Z])
        
        # Alternative: try direct point access for PolylineCurves
        elif hasattr(curve, 'PointCount'):
            point_count = curve.PointCount
            step = max(1, point_count // count)
            for i in range(0, point_count, step):
                if hasattr(curve, 'Point'):
                    pt = curve.Point(i)
                    if pt:
                        points.append([pt.X, pt.Y, pt.Z])
        
        # For PolylineCurves, try accessing points directly
        elif hasattr(curve, 'GetPoints'):
            curve_points = curve.GetPoints()
            if curve_points:
                step = max(1, len(curve_points) // count)
                for i in range(0, len(curve_points), step):
                    pt = curve_points[i]
                    points.append([pt.X, pt.Y, pt.Z])
                        
    except Exception as e:
        print(f"    ❌ Error sampling curve: {e}")
        
    return points

def extract_actual_geometry_from_3dm(model_file):
    """Extract real geometry from 3DM file - filtered by specific layers"""
    try:
        import rhino3dm as rhino
        print(f"📖 Reading 3DM file: {model_file}")
        
        model = rhino.File3dm.Read(model_file)
        if not model:
            print("❌ Failed to read 3DM file")
            return None
        
        print(f"📋 Found {len(model.Objects)} objects in the model")
        
        # Get layer information
        layers_by_index = {}
        if hasattr(model, 'Layers'):
            for layer in model.Layers:
                layers_by_index[layer.Index] = layer.Name
                print(f"  Layer {layer.Index}: '{layer.Name}'")
        
        # Target layers to extract
        target_layers = ["Costilla", "Vertebra", "Vias"]
        print(f"🎯 Target layers: {target_layers}")
        
        # Print object info for debugging with layer information
        relevant_objects = []
        for i, obj in enumerate(model.Objects):
            geom = obj.Geometry
            geom_type = type(geom).__name__ if geom else "None"
            
            # Get layer name
            layer_name = "Unknown"
            if hasattr(obj, 'Attributes') and hasattr(obj.Attributes, 'LayerIndex'):
                layer_index = obj.Attributes.LayerIndex
                layer_name = layers_by_index.get(layer_index, f"Layer_{layer_index}")
            
            print(f"  Object {i}: {geom_type} on layer '{layer_name}'")
            
            # Filter by target layers
            if layer_name in target_layers:
                relevant_objects.append((i, obj))
                print(f"    ✅ Will process this object")
        
        print(f"🔍 Found {len(relevant_objects)} objects on target layers")
        
        if len(relevant_objects) == 0:
            print("❌ No objects found on target layers")
            return None
        
        all_curves_with_layers = []
        
        # Process only objects on target layers
        for i, obj in relevant_objects:
            geom = obj.Geometry
            if not geom:
                continue
                
            # Get layer name for this object
            layer_name = "Unknown"
            if hasattr(obj, 'Attributes') and hasattr(obj.Attributes, 'LayerIndex'):
                layer_index = obj.Attributes.LayerIndex
                layer_name = layers_by_index.get(layer_index, f"Layer_{layer_index}")
                
            print(f"🔍 Processing object {i} ({type(geom).__name__}) on layer '{layer_name}'")
            
            try:
                # Handle different geometry types
                if hasattr(geom, 'Faces') and len(geom.Faces) > 0:
                    # It's a Brep/Surface - extract edge curves
                    print(f"  📐 Brep with {len(geom.Faces)} faces")
                    
                    # Try multiple approaches for Brep edge extraction
                    edges_extracted = 0
                    
                    # Method 1: Extract from Brep.Edges
                    if hasattr(geom, 'Edges'):
                        print(f"    🔍 Found {len(geom.Edges)} edges in Brep")
                        for edge_idx, edge in enumerate(geom.Edges):
                            try:
                                # Try different ways to get curve from edge
                                curve = None
                                
                                # Try various rhino3dm edge curve access methods
                                if hasattr(edge, 'EdgeCurve') and edge.EdgeCurve:
                                    curve = edge.EdgeCurve
                                elif hasattr(edge, 'Curve') and edge.Curve:
                                    curve = edge.Curve
                                elif hasattr(edge, 'DuplicateCurve'):
                                    curve = edge.DuplicateCurve()
                                elif hasattr(edge, 'ToNurbsCurve'):
                                    curve = edge.ToNurbsCurve()
                                
                                if curve:
                                    points = sample_curve_points(curve)
                                    if len(points) > 1:
                                        all_curves_with_layers.append((points, layer_name))
                                        edges_extracted += 1
                                        print(f"    ✅ Extracted edge {edge_idx} with {len(points)} points")
                                else:
                                    # Debug: print edge properties
                                    if edge_idx < 3:  # Only print for first few edges to avoid spam
                                        edge_props = [attr for attr in dir(edge) if not attr.startswith('_')]
                                        print(f"    🔍 Edge {edge_idx} properties: {edge_props[:10]}...")
                                    
                            except Exception as e:
                                print(f"    ⚠️ Edge {edge_idx} error: {e}")
                                if edge_idx < 3:  # Only print details for first few
                                    import traceback
                                    traceback.print_exc()
                    
                    # Method 2: Extract from Face loops if edges didn't work
                    if edges_extracted == 0:
                        print(f"    🔍 Trying face loop extraction...")
                        for face_idx, face in enumerate(geom.Faces):
                            try:
                                # Try different face boundary extraction methods
                                extracted_from_face = False
                                
                                # Method 2a: OuterLoop.To3dCurve()
                                if hasattr(face, 'OuterLoop'):
                                    loop = face.OuterLoop
                                    if hasattr(loop, 'To3dCurve'):
                                        curve = loop.To3dCurve()
                                        if curve:
                                            points = sample_curve_points(curve)
                                            if len(points) > 1:
                                                all_curves_with_layers.append((points, layer_name))
                                                edges_extracted += 1
                                                extracted_from_face = True
                                                print(f"    ✅ Extracted face {face_idx} outer loop with {len(points)} points")
                                
                                # Method 2b: Try getting surface boundary
                                if not extracted_from_face and hasattr(face, 'ToBrep'):
                                    try:
                                        face_brep = face.ToBrep()
                                        if face_brep and hasattr(face_brep, 'Edges'):
                                            for edge in face_brep.Edges:
                                                if hasattr(edge, 'DuplicateCurve'):
                                                    curve = edge.DuplicateCurve()
                                                    if curve:
                                                        points = sample_curve_points(curve)
                                                        if len(points) > 1:
                                                            all_curves_with_layers.append((points, layer_name))
                                                            edges_extracted += 1
                                                            extracted_from_face = True
                                                            break
                                    except:
                                        pass
                                
                                # Method 2c: Extract all loops
                                if not extracted_from_face and hasattr(face, 'Loops'):
                                    for loop_idx, loop in enumerate(face.Loops):
                                        if hasattr(loop, 'To3dCurve'):
                                            curve = loop.To3dCurve()
                                            if curve:
                                                points = sample_curve_points(curve)
                                                if len(points) > 1:
                                                    all_curves_with_layers.append((points, layer_name))
                                                    edges_extracted += 1
                                                    extracted_from_face = True
                                                    print(f"    ✅ Extracted face {face_idx} loop {loop_idx} with {len(points)} points")
                                
                                # Debug: print face properties if nothing worked
                                if not extracted_from_face and face_idx < 2:
                                    face_props = [attr for attr in dir(face) if not attr.startswith('_')]
                                    print(f"    🔍 Face {face_idx} properties: {face_props[:10]}...")
                                    
                            except Exception as e:
                                print(f"    ⚠️ Face {face_idx} error: {e}")
                                if face_idx < 2:
                                    import traceback
                                    traceback.print_exc()
                    
                    print(f"    📊 Total extracted from Brep: {edges_extracted} curves")
                
                elif hasattr(geom, 'ToBrep') or 'Extrusion' in type(geom).__name__:
                    # It's an Extrusion - extract profile curves using Profile3d method
                    print(f"  🏗️ Extrusion object - extracting profile curves")
                    edges_extracted = 0
                    try:
                        # Extract profile curves using Profile3d method
                        if hasattr(geom, 'Profile3d'):
                            try:
                                # Extract profile curves with different indices (start and end of extrusion)
                                for profile_idx in range(5):  # Try first few profiles
                                    try:
                                        # Get profile curve at parameter s=0.0 (start of extrusion)
                                        profile_curve = geom.Profile3d(profile_idx, 0.0)
                                        if profile_curve:
                                            points = sample_curve_points(profile_curve)
                                            if len(points) > 1:
                                                all_curves_with_layers.append((points, layer_name))
                                                edges_extracted += 1
                                        
                                        # Also try at parameter s=1.0 (end of extrusion) 
                                        profile_curve_end = geom.Profile3d(profile_idx, 1.0)
                                        if profile_curve_end:
                                            points_end = sample_curve_points(profile_curve_end)
                                            if len(points_end) > 1:
                                                all_curves_with_layers.append((points_end, layer_name))
                                                edges_extracted += 1
                                                
                                    except Exception:
                                        # Expected to fail when we run out of valid profiles
                                        break
                                        
                            except Exception as e:
                                print(f"    ❌ Profile3d access failed: {e}")
                        
                        # Try to extract path curves if available
                        for method in ['PathCurve', 'GetPath']:
                            if hasattr(geom, method):
                                try:
                                    value = getattr(geom, method)
                                    if callable(value):
                                        result = value()
                                        if result and not hasattr(result, 'X'):  # Not a point, try as curve
                                            points = sample_curve_points(result)
                                            if len(points) > 1:
                                                all_curves_with_layers.append((points, layer_name))
                                                edges_extracted += 1
                                except Exception:
                                    continue
                        
                        
                        # If direct access doesn't work, try converting to Brep as fallback
                        if edges_extracted == 0:
                            print(f"    🔄 Profile curves not found, trying Brep conversion")
                            brep = geom.ToBrep(True)
                            if brep and hasattr(brep, 'Edges'):
                                print(f"    🔍 Found {len(brep.Edges)} edges in Extrusion Brep")
                                
                                for edge_idx, edge in enumerate(brep.Edges):
                                    try:
                                        curve = None
                                        # Try alternative edge curve extraction methods
                                        if hasattr(edge, 'TryGetPlane'):
                                            # Skip planar edges for cleaner wireframes
                                            continue
                                        
                                        # Try direct curve access
                                        for method in ['DuplicateCurve', 'EdgeCurve', 'Curve']:
                                            if hasattr(edge, method):
                                                curve = getattr(edge, method)
                                                if callable(curve):
                                                    curve = curve()
                                                if curve:
                                                    break
                                        
                                        if curve:
                                            points = sample_curve_points(curve)
                                            if len(points) > 1:
                                                all_curves_with_layers.append((points, layer_name))
                                                edges_extracted += 1
                                                print(f"    ✅ Extracted edge {edge_idx} with {len(points)} points")
                                    except Exception as edge_e:
                                        continue
                        
                        print(f"    📊 Total extracted from Extrusion: {edges_extracted} curves")
                        
                    except Exception as extrusion_e:
                        print(f"  ❌ Error processing Extrusion: {extrusion_e}")
                        
                elif hasattr(geom, 'ToNurbsCurve') or 'Curve' in type(geom).__name__:
                    # It's a curve - sample it directly
                    points = sample_curve_points(geom)
                    if len(points) > 1:
                        all_curves_with_layers.append((points, layer_name))
                        print(f"  ✅ Extracted curve with {len(points)} points")
                
                elif hasattr(geom, 'Vertices'):
                    # It's a mesh
                    print(f"  📐 Mesh with {len(geom.Vertices)} vertices")
                    # Extract mesh edges as wireframe
                    mesh_curves = extract_mesh_edges(geom)
                    # Add layer info to mesh curves
                    for curve in mesh_curves:
                        all_curves_with_layers.append((curve, layer_name))
                    print(f"  ✅ Extracted {len(mesh_curves)} mesh edge curves")
                
            except Exception as e:
                print(f"  ❌ Error processing object {i}: {e}")
                continue
        
        if len(all_curves_with_layers) > 0:
            print(f"🎉 Successfully extracted {len(all_curves_with_layers)} curves from target layers")
            return all_curves_with_layers
        else:
            print("❌ No curves could be extracted from target layers")
            print("💡 This could mean:")
            print("   - Layer names don't match exactly")
            print("   - Objects don't have extractable geometry")
            print("   - API calls are failing")
            return None
            
    except ImportError:
        print("❌ rhino3dm library not available")
        return None
    except Exception as e:
        print("❌ Error reading 3DM file: {e}")
        import traceback
        traceback.print_exc()
        return None

def extract_mesh_edges(mesh):
    """Extract edge curves from mesh"""
    curves = []
    try:
        if hasattr(mesh, 'Faces') and hasattr(mesh, 'Vertices'):
            # Simple approach: extract face boundaries
            for face in mesh.Faces:
                if hasattr(face, 'IsQuad') and face.IsQuad:
                    # Quad face
                    indices = [face.A, face.B, face.C, face.D, face.A]  # Close the loop
                else:
                    # Triangle face
                    indices = [face.A, face.B, face.C, face.A]  # Close the loop
                
                points = []
                for idx in indices:
                    if idx < len(mesh.Vertices):
                        v = mesh.Vertices[idx]
                        points.append([v.X, v.Y, v.Z])
                
                if len(points) > 1:
                    curves.append(points)
                    
    except Exception as e:
        print(f"    ❌ Error extracting mesh edges: {e}")
        
    return curves

def calculate_bounds(wireframe_lines_with_layers):
    """Calculate bounding box from wireframe curves with layer info"""
    all_points = []
    for line_data in wireframe_lines_with_layers:
        if isinstance(line_data, tuple):
            line, layer_name = line_data
        else:
            # Backward compatibility
            line = line_data
        all_points.extend(line)
    
    if not all_points:
        return (-35, -35, -12, 35, 35, 12)  # Default bounds
    
    points_array = np.array(all_points)
    padding = 5  # Add some padding
    
    return (
        points_array[:, 0].min() - padding, 
        points_array[:, 1].min() - padding, 
        points_array[:, 2].min() - padding,
        points_array[:, 0].max() + padding, 
        points_array[:, 1].max() + padding, 
        points_array[:, 2].max() + padding
    )

def create_wireframe_view(wireframe_lines_with_layers, view_name, bounds, output_dir='images'):
    """Create CAD-style wireframe view with layer colors"""
    print(f"🖼️  Creating {view_name} view with {len(wireframe_lines_with_layers)} curves...")
    
    min_x, min_y, min_z, max_x, max_y, max_z = bounds
    
    # Define colors for each layer
    layer_colors = {
        'Costilla': '#2E86AB',    # Professional blue
        'Vertebra': '#A23B72',    # Deep magenta/red
        'Vias': '#F18F01'        # Professional orange
    }
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7.5), dpi=150)
    ax.set_facecolor('white')
    
    # Plot each wireframe curve with layer-specific color
    for line_data in wireframe_lines_with_layers:
        if isinstance(line_data, tuple):
            line, layer_name = line_data
            color = layer_colors.get(layer_name, '#333333')  # Default gray
        else:
            # Backward compatibility for old format
            line = line_data
            color = '#333333'
            
        if len(line) < 2:
            continue
            
        line_array = np.array(line)
        
        # Project onto appropriate plane
        if view_name == 'top':
            x, y = line_array[:, 0], line_array[:, 1]
            ax.plot(x, y, color=color, linewidth=1.2, alpha=0.85)
            
        elif view_name == 'front':
            x, y = line_array[:, 0], line_array[:, 2]
            ax.plot(x, y, color=color, linewidth=1.2, alpha=0.85)
            
        elif view_name == 'right':
            x, y = line_array[:, 1], line_array[:, 2]
            ax.plot(x, y, color=color, linewidth=1.2, alpha=0.85)
            
        elif view_name == 'isometric':
            # Isometric projection with clean wireframes and depth styling
            iso_x = line_array[:, 0] - line_array[:, 1] * 0.5
            iso_y = line_array[:, 2] + (line_array[:, 0] + line_array[:, 1]) * 0.25

            # Calculate depth for clean wireframe styling
            depth = line_array[:, 2] + line_array[:, 1] * 0.3
            depth_normalized = (depth - np.min(depth)) / (np.max(depth) - np.min(depth) + 1e-10)

            # Varying line width and alpha based on depth for clean 3D effect
            base_alpha = 0.9
            depth_alpha = base_alpha - depth_normalized * 0.3  # Front lines more opaque
            avg_alpha = np.mean(depth_alpha)

            base_width = 1.0
            depth_width = base_width + (1 - depth_normalized) * 0.5  # Front lines thicker
            avg_width = np.mean(depth_width)

            ax.plot(iso_x, iso_y, color=color, linewidth=avg_width, alpha=avg_alpha)
    
    # Technical CAD-style wireframe views with layer colors
    
    # Set view-specific properties with grids and axes (no titles)
    if view_name == 'top':
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
        ax.set_xlabel('X', fontsize=12, fontweight='bold')
        ax.set_ylabel('Y', fontsize=12, fontweight='bold')
        
    elif view_name == 'front':
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_z, max_z)
        ax.set_xlabel('X', fontsize=12, fontweight='bold')
        ax.set_ylabel('Z', fontsize=12, fontweight='bold')
        
    elif view_name == 'right':
        ax.set_xlim(min_y, max_y)
        ax.set_ylim(min_z, max_z)
        ax.set_xlabel('Y', fontsize=12, fontweight='bold')
        ax.set_ylabel('Z', fontsize=12, fontweight='bold')
        
    elif view_name == 'isometric':
        # Calculate bounds for isometric view with more padding (zoom out)
        iso_bounds_x = []
        iso_bounds_y = []
        for line_data in wireframe_lines_with_layers:
            if isinstance(line_data, tuple):
                line, layer_name = line_data
            else:
                line = line_data
            if len(line) >= 2:
                line_array = np.array(line)
                iso_x = line_array[:, 0] - line_array[:, 1] * 0.5
                iso_y = line_array[:, 2] + (line_array[:, 0] + line_array[:, 1]) * 0.25
                iso_bounds_x.extend(iso_x)
                iso_bounds_y.extend(iso_y)
        
        if iso_bounds_x and iso_bounds_y:
            padding = 25  # Increased padding for zoom out effect
            ax.set_xlim(min(iso_bounds_x) - padding, max(iso_bounds_x) + padding)
            ax.set_ylim(min(iso_bounds_y) - padding, max(iso_bounds_y) + padding)
        
        ax.set_xlabel('X-Y', fontsize=12, fontweight='bold')
        ax.set_ylabel('Z+(X+Y)/4', fontsize=12, fontweight='bold')
    
    # Professional CAD-style appearance
    ax.set_aspect('equal', adjustable='box')
    
    # Add grid and centerlines for technical drawings
    ax.grid(True, alpha=0.3, linewidth=0.5, color='gray')
    
    # Add centerlines for all views except isometric
    if view_name != 'isometric':
        ax.axhline(0, color='red', linestyle='--', alpha=0.4, linewidth=1)
        ax.axvline(0, color='red', linestyle='--', alpha=0.4, linewidth=1)
    
    # Professional axis styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_color('gray')
    
    # Save
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f'{view_name}.png')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', 
               facecolor='white', edgecolor='none', format='png')
    plt.close()
    print(f"✅ Saved {view_name} view to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Generate 3D views from Rhino 3DM files')
    parser.add_argument('--input', '-i', help='Input 3DM file path')
    parser.add_argument('--output', '-o', default='images', help='Output directory for images')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    print("🚀 Starting 3D view generation...")
    
    # Find 3DM file
    if args.input:
        model_file = args.input
        if not os.path.exists(model_file):
            print(f"❌ Input file not found: {model_file}")
            return False
    else:
        model_file = find_3dm_file()
        if not model_file:
            print("❌ No .3dm file found in current directory")
            return False
    
    print(f"📁 Using 3DM file: {model_file}")
    
    # Extract geometry from specific layers only
    wireframe_lines = extract_actual_geometry_from_3dm(model_file)
    
    # Fail if no geometry was extracted
    if wireframe_lines is None or len(wireframe_lines) == 0:
        print("💥 FAILED: No geometry could be extracted from the 3DM file")
        print("🔍 Check:")
        print("   1. Layer names match exactly: 'Costilla' 'Vertebra' and 'Vias'")
        print("   2. Objects on those layers have valid geometry")
        print("   3. rhino3dm library can access the geometry")
        return False
    
    # Calculate bounds
    bounds = calculate_bounds(wireframe_lines)
    print(f"📏 Geometry bounds: {bounds}")
    
    # Generate views
    views = ['top', 'front', 'right', 'isometric']
    for view in views:
        create_wireframe_view(wireframe_lines, view, bounds, args.output)
    
    print("🎉 All views generated successfully from real geometry!")
    return True

if __name__ == "__main__":
    try:
        success = main()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"💥 Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1) 