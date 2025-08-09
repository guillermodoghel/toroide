#!/bin/bash
# Setup script for local development

echo "🚀 Setting up 3D model view generation environment..."

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.8+ first."
    exit 1
fi

# Create virtual environment
echo "📦 Creating virtual environment..."
python3 -m venv venv

# Activate virtual environment
echo "🔌 Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "⬆️  Upgrading pip..."
pip install --upgrade pip

# Install dependencies
echo "📋 Installing dependencies from requirements.txt..."
pip install -r requirements.txt

echo "✅ Setup complete!"
echo ""
echo "To use the environment:"
echo "  source venv/bin/activate"
echo ""
echo "To generate views:"
echo "  python generate_3d_views.py --output images"
echo ""
echo "To test with fallback torus:"
echo "  python generate_3d_views.py --fallback"
echo ""
echo "To deactivate when done:"
echo "  deactivate" 