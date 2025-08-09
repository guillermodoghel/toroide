# Toroide 3D Model

This repository contains a 3D model of a toroide (torus/donut shape) created in Rhino 3D.

## 3D Model Views

<div align="center">

### Isometric View
![Isometric View](images/isometric.png)

<table>
  <tr>
    <td align="center">
      <strong>Top View</strong><br/>
      <img src="images/top.png" alt="Top View" width="300"/>
    </td>
    <td align="center">
      <strong>Front View</strong><br/>
      <img src="images/front.png" alt="Front View" width="300"/>
    </td>
    <td align="center">
      <strong>Right View</strong><br/>
      <img src="images/right.png" alt="Right View" width="300"/>
    </td>
  </tr>
</table>

</div>

## Interactive 3D Viewing

To view this model interactively:

1. **Download the file**: [`${MODEL_NAME}`](${MODEL_NAME}) (${MODEL_SIZE})
2. **Online viewers**:
   - [3DViewer.net](https://3dviewer.net/) - Upload and view in browser
   - [Online 3D Viewer](https://viewer.3dprintcloud.com/) - Another web-based option
3. **Desktop software**:
   - Rhino 3D (native format)
   - FreeCAD (open source)
   - Blender (with import plugins)

## Model Information

| Property | Value |
|----------|-------|
| **Format** | Rhino 3D (.3dm) |
| **File Size** | ${MODEL_SIZE} |
| **Last Updated** | ${CURRENT_DATE} |
| **Views Generated** | ${CURRENT_DATETIME} |
| **Type** | Parametric torus/donut shape |

## Files Structure

```
toroide/
├── ${MODEL_NAME}           # Main 3D model file
├── images/                      # Rendered view images
│   ├── front.png               # Front orthographic view
│   ├── isometric.png           # Isometric (3/4) view
│   ├── right.png               # Right side orthographic view
│   └── top.png                 # Top orthographic view
└── README.md                   # This documentation
```

## Technical Notes

- The model is created and maintained in **Rhino 3D**
- Views are automatically generated using **rhino3dm** Python library
- Static orthographic and isometric views are rendered as PNG images
- For the best interactive experience, download the `.3dm` file and open in compatible software
- GitHub doesn't support native 3D file preview, but the static images provide comprehensive views

## Generation Details

- **View Generation**: Automated via GitHub Actions using rhino3dm library
- **Rendering**: matplotlib with custom orthographic projections
- **Update Frequency**: Automatically triggered on model file changes
- **Fallback**: Creates sample torus if original geometry cannot be processed

---

<div align="center">
<em>Documentation and views automatically updated on ${CURRENT_DATETIME}</em><br/>
<em>🔄 Generated using rhino3dm Python library</em>
</div>
