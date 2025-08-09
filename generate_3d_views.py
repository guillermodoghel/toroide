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

def extract_actual_geometry_from_3dm(model_file):
    """Extract real geometry from 3DM file - not just bounding boxes"""
    try:
        import rhino3dm as rhino
        print(f"📖 Reading 3DM file: {model_file}")
        
        model = rhino.File3dm.Read(model_file)
        if not model:
            print("❌ Failed to read 3DM file")
            return None
        
        print(f"📋 Found {len(model.Objects)} objects in the model")
        
        # Print object info for debugging
        for i, obj in enumerate(model.Objects):
            geom = obj.Geometry
            geom_type = type(geom).__name__ if geom else "None"
            print(f"  Object {i}: {geom_type}")
        
        all_curves = []
        
        # Process each object
        for i, obj in enumerate(model.Objects):
            geom = obj.Geometry
            if not geom:
                continue
                
            print(f"🔍 Processing object {i} ({type(geom).__name__})")
            
            try:
                # Handle different geometry types
                if hasattr(geom, 'Faces') and len(geom.Faces) > 0:
                    # It's a Brep/Surface - extract edge curves
                    print(f"  📐 Brep with {len(geom.Faces)} faces")
                    for face in geom.Faces:
                        # Get outer boundary
                        if hasattr(face, 'OuterLoop'):
                            loop = face.OuterLoop
                            if hasattr(loop, 'To3dCurve'):
                                curve = loop.To3dCurve()
                                if curve:
                                    points = sample_curve_points(curve)
                                    if len(points) > 1:
                                        all_curves.append(points)
                                        print(f"    ✅ Extracted boundary curve with {len(points)} points")
                        
                        # Get inner holes/loops
                        if hasattr(face, 'Loops'):
                            for loop in face.Loops:
                                if hasattr(loop, 'To3dCurve'):
                                    curve = loop.To3dCurve()
                                    if curve:
                                        points = sample_curve_points(curve)
                                        if len(points) > 1:
                                            all_curves.append(points)
                                            print(f"    ✅ Extracted loop curve with {len(points)} points")
                
                elif hasattr(geom, 'ToNurbsCurve'):
                    # It's a curve
                    curve = geom.ToNurbsCurve() if callable(geom.ToNurbsCurve) else geom.ToNurbsCurve
                    if curve:
                        points = sample_curve_points(curve)
                        if len(points) > 1:
                            all_curves.append(points)
                            print(f"  ✅ Extracted curve with {len(points)} points")
                
                elif hasattr(geom, 'Vertices'):
                    # It's a mesh
                    print(f"  📐 Mesh with {len(geom.Vertices)} vertices")
                    # Extract mesh edges as wireframe
                    mesh_curves = extract_mesh_edges(geom)
                    all_curves.extend(mesh_curves)
                    print(f"  ✅ Extracted {len(mesh_curves)} mesh edge curves")
                
                else:
                    # Try to convert to Brep first
                    if hasattr(geom, 'ToBrep'):
                        brep = geom.ToBrep()
                        if brep and hasattr(brep, 'Faces'):
                            print(f"  📐 Converted to Brep with {len(brep.Faces)} faces")
                            for face in brep.Faces:
                                if hasattr(face, 'OuterLoop'):
                                    loop = face.OuterLoop
                                    if hasattr(loop, 'To3dCurve'):
                                        curve = loop.To3dCurve()
                                        if curve:
                                            points = sample_curve_points(curve)
                                            if len(points) > 1:
                                                all_curves.append(points)
                                                print(f"    ✅ Extracted boundary curve with {len(points)} points")
                
            except Exception as e:
                print(f"  ❌ Error processing object {i}: {e}")
                continue
        
        if len(all_curves) > 0:
            print(f"🎉 Successfully extracted {len(all_curves)} curves from 3DM file")
            return all_curves
        else:
            print("❌ No curves could be extracted from 3DM file")
            return None
            
    except ImportError:
        print("❌ rhino3dm library not available")
        return None
    except Exception as e:
        print(f"❌ Error reading 3DM file: {e}")
        import traceback
        traceback.print_exc()
        return None

def sample_curve_points(curve, count=50):
    """Sample points along a curve"""
    points = []
    try:
        if hasattr(curve, 'Domain'):
            domain = curve.Domain
            for i in range(count):
                t = domain.Min + (domain.Max - domain.Min) * i / (count - 1)
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
                
    except Exception as e:
        print(f"    ❌ Error sampling curve: {e}")
        
    return points

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

def create_fallback_torus_wireframe(major_radius=25, minor_radius=8, u_count=32, v_count=16):
    """Create wireframe lines for a torus as fallback"""
    print("🔄 Creating fallback torus wireframe...")
    lines = []
    
    # U-direction curves (around the major radius)
    for i in range(u_count):
        u = 2 * np.pi * i / u_count
        curve_points = []
        for j in range(v_count + 1):  # +1 to close the curve
            v = 2 * np.pi * j / v_count
            x = (major_radius + minor_radius * np.cos(v)) * np.cos(u)
            y = (major_radius + minor_radius * np.cos(v)) * np.sin(u)
            z = minor_radius * np.sin(v)
            curve_points.append([x, y, z])
        lines.append(curve_points)
    
    # V-direction curves (around the minor radius)
    for j in range(0, v_count, 2):  # Every other line to avoid clutter
        v = 2 * np.pi * j / v_count
        curve_points = []
        for i in range(u_count + 1):  # +1 to close the curve
            u = 2 * np.pi * i / u_count
            x = (major_radius + minor_radius * np.cos(v)) * np.cos(u)
            y = (major_radius + minor_radius * np.cos(v)) * np.sin(u)
            z = minor_radius * np.sin(v)
            curve_points.append([x, y, z])
        lines.append(curve_points)
    
    return lines

def calculate_bounds(wireframe_lines):
    """Calculate bounding box from wireframe curves"""
    all_points = []
    for line in wireframe_lines:
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

def create_wireframe_view(wireframe_lines, view_name, bounds, output_dir='images'):
    """Create CAD-style wireframe view"""
    print(f"🖼️  Creating {view_name} view with {len(wireframe_lines)} curves...")
    
    min_x, min_y, min_z, max_x, max_y, max_z = bounds
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7.5), dpi=150)
    ax.set_facecolor('white')
    
    # Plot each wireframe curve
    for line in wireframe_lines:
        if len(line) < 2:
            continue
            
        line_array = np.array(line)
        
        # Project onto appropriate plane
        if view_name == 'top':
            x, y = line_array[:, 0], line_array[:, 1]
            ax.plot(x, y, 'k-', linewidth=1.0, alpha=0.8)
            
        elif view_name == 'front':
            x, y = line_array[:, 0], line_array[:, 2]
            ax.plot(x, y, 'k-', linewidth=1.0, alpha=0.8)
            
        elif view_name == 'right':
            x, y = line_array[:, 1], line_array[:, 2]
            ax.plot(x, y, 'k-', linewidth=1.0, alpha=0.8)
            
        elif view_name == 'isometric':
            # Isometric projection
            iso_x = line_array[:, 0] - line_array[:, 1] * 0.5
            iso_y = line_array[:, 2] + (line_array[:, 0] + line_array[:, 1]) * 0.25
            ax.plot(iso_x, iso_y, 'k-', linewidth=0.8, alpha=0.7)
    
    # Set view-specific properties
    if view_name == 'top':
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
        ax.set_xlabel('X', fontsize=12, fontweight='bold')
        ax.set_ylabel('Y', fontsize=12, fontweight='bold')
        ax.set_title('Top View', fontsize=16, pad=20, fontweight='bold')
        
    elif view_name == 'front':
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_z, max_z)
        ax.set_xlabel('X', fontsize=12, fontweight='bold')
        ax.set_ylabel('Z', fontsize=12, fontweight='bold')
        ax.set_title('Front View', fontsize=16, pad=20, fontweight='bold')
        
    elif view_name == 'right':
        ax.set_xlim(min_y, max_y)
        ax.set_ylim(min_z, max_z)
        ax.set_xlabel('Y', fontsize=12, fontweight='bold')
        ax.set_ylabel('Z', fontsize=12, fontweight='bold')
        ax.set_title('Right View', fontsize=16, pad=20, fontweight='bold')
        
    elif view_name == 'isometric':
        ax.set_title('Isometric View', fontsize=16, pad=20, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
    
    # Styling
    ax.set_aspect('equal', adjustable='box')
    
    if view_name != 'isometric':
        ax.grid(True, alpha=0.3, linewidth=0.5, color='gray')
        # Add centerlines
        ax.axhline(0, color='red', linestyle='--', alpha=0.4, linewidth=1)
        ax.axvline(0, color='red', linestyle='--', alpha=0.4, linewidth=1)
    
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
    parser.add_argument('--fallback', action='store_true', help='Force use fallback torus')
    
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
    
    # Extract geometry
    wireframe_lines = None
    if not args.fallback:
        wireframe_lines = extract_actual_geometry_from_3dm(model_file)
    
    # Use fallback if extraction failed
    if wireframe_lines is None or len(wireframe_lines) == 0:
        print("⚠️  Using fallback torus wireframe")
        wireframe_lines = create_fallback_torus_wireframe()
    
    # Calculate bounds
    bounds = calculate_bounds(wireframe_lines)
    print(f"📏 Geometry bounds: {bounds}")
    
    # Generate views
    views = ['top', 'front', 'right', 'isometric']
    for view in views:
        create_wireframe_view(wireframe_lines, view, bounds, args.output)
    
    print("🎉 All views generated successfully!")
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