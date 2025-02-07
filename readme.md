# 🧬 Protein Sequence Analysis Tool

A Python-based tool for analyzing protein sequences and calculating various physicochemical properties. This tool uses BioPython's ProtParam module to analyze FASTA sequences and exports the results to an Excel file.

## 🎯 Features

The tool calculates the following protein properties:
- Molecular Weight
- Amino Acid Composition
- Molar Extinction Coefficient
- Isoelectric Point
- Instability Index
- Aromaticity
- GRAVY (Grand Average of Hydropathy)
- Flexibility

## 🛠️ Prerequisites

- Python 3.x
- BioPython
- Pandas
- OpenPyXL

## ⚙️ Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/protein-sequence-analysis.git
cd protein-sequence-analysis
```

2. Install required packages:
```bash
pip install biopython pandas openpyxl
```

## 📋 Input Requirements

- Input file should be in FASTA format
- File should be named "Sequence.fasta"
- Protein sequence should contain valid amino acid letters

Example FASTA format:
```
>Protein_Name
MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQ
```

## 🚀 Usage

1. Place your FASTA file (named "Sequence.fasta") in the same directory as the script

2. Run the script:
```bash
python protein_analysis.py
```

3. Check the output file "protein_feature_data.xlsx" for results

## 📊 Output

The script generates an Excel file ("protein_feature_data.xlsx") with the following columns:
- Molecular_Weight
- Amino_Acid_Count
- molar_extinction_coefficient
- isoelectric_point
- instability_index
- aromaticity
- Gravy
- Flexibility

## 📝 Code Description

```python
# Read FASTA file
input_file = open("Sequence.fasta", "r")
for record in SeqIO.parse(input_file, "fasta"):
    # Process sequence and calculate properties
    my_sec = str(record.seq).rstrip('\\')
    analyse = ProteinAnalysis(my_sec)
    
    # Calculate various properties
    mol_weight = analyse.molecular_weight()
    count_amino = analyse.count_amino_acids()
    # ... other calculations
```

## 📊 Visualizations

You can enhance your analysis by visualizing the protein properties. Here are some example visualization codes using matplotlib and seaborn:

### 1. Amino Acid Composition Bar Plot
```python
import matplotlib.pyplot as plt
import seaborn as sns

def plot_amino_acid_composition(count_amino):
    plt.figure(figsize=(12, 6))
    sns.barplot(x=list(count_amino.keys()), y=list(count_amino.values()))
    plt.title('Amino Acid Composition')
    plt.xlabel('Amino Acids')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('amino_acid_composition.png')
    plt.close()

# Usage
plot_amino_acid_composition(count_amino)
```
![Amino Acid Composition Example](https://raw.githubusercontent.com/yourusername/protein-analysis/main/examples/amino_acid_plot.png)

### 2. Protein Flexibility Profile
```python
def plot_flexibility_profile(flex):
    plt.figure(figsize=(12, 6))
    plt.plot(range(len(flex)), flex, '-b')
    plt.title('Protein Flexibility Profile')
    plt.xlabel('Position')
    plt.ylabel('Flexibility Score')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('flexibility_profile.png')
    plt.close()

# Usage
plot_flexibility_profile(flex)
```
![Flexibility Profile Example](https://raw.githubusercontent.com/yourusername/protein-analysis/main/examples/flexibility_plot.png)

### 3. Property Comparison Radar Chart
```python
import numpy as np

def plot_property_radar(properties):
    # Normalize the properties
    normalized_props = {
        'MW': properties['Molecular_Weight'] / 100000,
        'pI': properties['isoelectric_point'] / 14,
        'Instability': properties['instability_index'] / 100,
        'Aromaticity': properties['aromaticity'],
        'GRAVY': (properties['Gravy'] + 4) / 8  # Normalize between -4 and 4
    }
    
    categories = list(normalized_props.keys())
    values = list(normalized_props.values())
    
    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False)
    values = np.concatenate((values, [values[0]]))
    angles = np.concatenate((angles, [angles[0]]))
    
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    ax.plot(angles, values)
    ax.fill(angles, values, alpha=0.25)
    ax.set_thetagrids(angles[:-1] * 180/np.pi, categories)
    plt.title('Protein Properties Radar Chart')
    plt.tight_layout()
    plt.savefig('property_radar.png')
    plt.close()

# Usage
properties = {
    'Molecular_Weight': mol_weight,
    'isoelectric_point': iso_point,
    'instability_index': ist_index,
    'aromaticity': aromati,
    'Gravy': gra_vy
}
plot_property_radar(properties)
```
![Property Radar Example](https://raw.githubusercontent.com/yourusername/protein-analysis/main/examples/radar_plot.png)

### 4. Combined Analysis Dashboard
```python
def create_analysis_dashboard(count_amino, flex, properties):
    fig = plt.figure(figsize=(15, 10))
    
    # Amino Acid Composition
    ax1 = plt.subplot(221)
    sns.barplot(x=list(count_amino.keys()), y=list(count_amino.values()), ax=ax1)
    ax1.set_title('Amino Acid Composition')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
    
    # Flexibility Profile
    ax2 = plt.subplot(222)
    ax2.plot(range(len(flex)), flex, '-b')
    ax2.set_title('Flexibility Profile')
    
    # Property Radar
    ax3 = plt.subplot(223, projection='polar')
    categories = list(properties.keys())
    values = list(properties.values())
    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False)
    ax3.plot(angles, values)
    ax3.fill(angles, values, alpha=0.25)
    ax3.set_title('Property Radar')
    
    # Save dashboard
    plt.tight_layout()
    plt.savefig('protein_analysis_dashboard.png')
    plt.close()

# Usage
create_analysis_dashboard(count_amino, flex, properties)
```
![Analysis Dashboard Example](https://raw.githubusercontent.com/yourusername/protein-analysis/main/examples/dashboard.png)

### Additional Visualization Tips:

1. **Data Export for External Tools**:
```python
# Export data for external visualization
def export_visualization_data(count_amino, flex, properties):
    # Create separate CSV files for each analysis
    pd.DataFrame(count_amino.items(), 
                columns=['Amino_Acid', 'Count']).to_csv('amino_acid_data.csv')
    pd.DataFrame({'Position': range(len(flex)), 
                 'Flexibility': flex}).to_csv('flexibility_data.csv')
    pd.DataFrame(properties, index=[0]).to_csv('properties_data.csv')
```

2. **Interactive Visualization with Plotly**:
```python
import plotly.express as px
import plotly.graph_objects as go

def create_interactive_plot(count_amino):
    fig = px.bar(x=list(count_amino.keys()), 
                 y=list(count_amino.values()),
                 title='Interactive Amino Acid Composition')
    fig.write_html('interactive_plot.html')
```

## 📊 Required Visualization Packages

Add these to your requirements.txt:
```
matplotlib>=3.5.0
seaborn>=0.11.0
plotly>=5.3.0  # For interactive plots
```

## 📌 Important Notes

- Ensure your FASTA sequence contains only valid amino acid letters
- The tool processes one sequence at a time
- Large sequences may take longer to process
- Output will overwrite existing Excel file with the same name

## 🐛 Troubleshooting

Common issues and solutions:

1. Invalid Sequence Error:
   - Ensure sequence contains only valid amino acid letters
   - Remove any special characters or spaces

2. File Not Found Error:
   - Check if "Sequence.fasta" exists in the correct directory
   - Verify file name and extension

3. Excel File Access Error:
   - Close the output Excel file before running the script
   - Check write permissions in the directory

## 🤝 Contributing

1. Fork the repository
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Create a new Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Contact

For questions or support, please open an issue in the GitHub repository.

---
*Made with ❤️ for the bioinformatics community*
