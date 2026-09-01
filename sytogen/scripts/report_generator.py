"""
Report generation for SyToGen results and DNA analysis.

Generates:
- HTML reports with interactive visualizations
- JSON summaries for programmatic access
- Markdown reports for documentation
"""

import json
import base64
from datetime import datetime
from typing import Dict, List, Optional, Any
from pathlib import Path


class ReportBuilder:
    """Build comprehensive reports for SyToGen analysis."""
    
    def __init__(self, title: str = "SyToGen Report"):
        self.title = title
        self.sections = []
        self.metadata = {
            "generated": datetime.now().isoformat(),
            "title": title,
        }
    
    def add_section(self, title: str, content: str, section_type: str = "text"):
        """Add a section to the report."""
        self.sections.append({
            "title": title,
            "content": content,
            "type": section_type,
        })
    
    def add_metric(self, label: str, value: Any, unit: str = ""):
        """Add a key metric to the report."""
        if not hasattr(self, "metrics"):
            self.metrics = []
        self.metrics.append({
            "label": label,
            "value": value,
            "unit": unit,
        })
    
    def add_chart_data(self, name: str, data: Dict[str, Any]):
        """Add chart-ready data."""
        if not hasattr(self, "charts"):
            self.charts = {}
        self.charts[name] = data
    
    def to_html(self, include_css: bool = True) -> str:
        """Generate HTML report."""
        html = []
        
        if include_css:
            html.append(self._get_css())
        
        html.append(f"""
        <html>
        <head>
            <meta charset="UTF-8">
            <title>{self.title}</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        </head>
        <body>
            <div class="report-container">
                <h1>{self.title}</h1>
                <div class="metadata">
                    <p>Generated: {self.metadata['generated']}</p>
                </div>
        """)
        
        # Metrics dashboard
        if hasattr(self, "metrics") and self.metrics:
            html.append('<div class="metrics-grid">')
            for metric in self.metrics:
                val_str = f"{metric['value']:.1f}" if isinstance(metric['value'], float) else str(metric['value'])
                html.append(f"""
                <div class="metric-card">
                    <div class="metric-label">{metric['label']}</div>
                    <div class="metric-value">{val_str} {metric['unit']}</div>
                </div>
                """)
            html.append('</div>')
        
        # Charts
        if hasattr(self, "charts"):
            for chart_name, chart_data in self.charts.items():
                html.append(f'<div id="chart-{chart_name}"></div>')
                html.append(self._chart_js(chart_name, chart_data))
        
        # Sections
        for section in self.sections:
            html.append(f'<div class="section">')
            html.append(f'<h2>{section["title"]}</h2>')
            if section["type"] == "table":
                html.append(f'<div class="table-container">{section["content"]}</div>')
            elif section["type"] == "code":
                html.append(f'<pre><code>{self._escape_html(section["content"])}</code></pre>')
            else:
                html.append(f'<div class="section-content">{self._escape_html(section["content"])}</div>')
            html.append('</div>')
        
        html.append("""
            </div>
        </body>
        </html>
        """)
        
        return "\n".join(html)
    
    def to_json(self) -> str:
        """Generate JSON report."""
        report_data = {
            "metadata": self.metadata,
            "sections": self.sections,
        }
        
        if hasattr(self, "metrics"):
            report_data["metrics"] = self.metrics
        
        if hasattr(self, "charts"):
            report_data["charts"] = self.charts
        
        return json.dumps(report_data, indent=2)
    
    def to_markdown(self) -> str:
        """Generate Markdown report."""
        md = []
        md.append(f"# {self.title}\n")
        md.append(f"**Generated**: {self.metadata['generated']}\n\n")
        
        if hasattr(self, "metrics"):
            md.append("## Key Metrics\n")
            for metric in self.metrics:
                val = f"{metric['value']:.2f}" if isinstance(metric['value'], float) else metric['value']
                md.append(f"- **{metric['label']}**: {val} {metric['unit']}\n")
            md.append("\n")
        
        for section in self.sections:
            md.append(f"## {section['title']}\n")
            md.append(f"{section['content']}\n\n")
        
        return "".join(md)
    
    def _get_css(self) -> str:
        """Return embedded CSS for HTML reports."""
        return """
        <style>
            body {
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background: #f5f5f5;
                margin: 0;
                padding: 20px;
            }
            .report-container {
                max-width: 1200px;
                margin: 0 auto;
                background: white;
                padding: 40px;
                border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }
            h1 {
                color: #2c3e50;
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
            }
            h2 {
                color: #34495e;
                margin-top: 30px;
            }
            .metadata {
                font-size: 0.9em;
                color: #7f8c8d;
                font-style: italic;
            }
            .metrics-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
                margin: 20px 0;
            }
            .metric-card {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            }
            .metric-label {
                font-size: 0.9em;
                opacity: 0.9;
            }
            .metric-value {
                font-size: 2em;
                font-weight: bold;
                margin-top: 10px;
            }
            .section {
                margin: 30px 0;
                padding: 20px;
                background: #f9f9f9;
                border-left: 4px solid #3498db;
                border-radius: 4px;
            }
            .section-content {
                line-height: 1.6;
                color: #333;
            }
            .table-container {
                overflow-x: auto;
            }
            table {
                width: 100%;
                border-collapse: collapse;
                margin: 10px 0;
            }
            th, td {
                padding: 10px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }
            th {
                background: #34495e;
                color: white;
            }
            tr:hover {
                background: #f0f0f0;
            }
            pre {
                background: #2c3e50;
                color: #ecf0f1;
                padding: 15px;
                border-radius: 4px;
                overflow-x: auto;
            }
            code {
                font-family: 'Courier New', monospace;
            }
        </style>
        """
    
    def _chart_js(self, name: str, data: Dict) -> str:
        """Generate Plotly.js chart initialization."""
        data_json = json.dumps(data.get("data", []))
        layout_json = json.dumps(data.get("layout", {}))
        
        return f"""
        <script>
            var data_{name} = {data_json};
            var layout_{name} = {layout_json};
            Plotly.newPlot('chart-{name}', data_{name}, layout_{name}, {{responsive: true}});
        </script>
        """
    
    @staticmethod
    def _escape_html(text: str) -> str:
        """Escape HTML special characters."""
        return (text.replace("&", "&amp;")
                   .replace("<", "&lt;")
                   .replace(">", "&gt;")
                   .replace('"', "&quot;"))


def create_sequence_analysis_report(sequence: str, quality_scores: Dict[str, int],
                                   repeats: List, structures: List) -> ReportBuilder:
    """
    Create a report for sequence quality analysis.
    
    Args:
        sequence: DNA sequence analyzed
        quality_scores: Dict of quality metrics
        repeats: List of detected repeats
        structures: List of detected secondary structures
        
    Returns:
        ReportBuilder with formatted analysis
    """
    report = ReportBuilder(f"Sequence Analysis Report - {len(sequence)}bp")
    
    # Add quality metrics
    for metric_name, score in quality_scores.items():
        if metric_name != "overall":
            report.add_metric(
                metric_name.replace("_", " ").title(),
                score,
                "%"
            )
    
    report.add_metric("Overall Quality", quality_scores.get("overall", 0), "%")
    
    # Quality summary
    report.add_section(
        "Summary",
        f"""
        <p>Analyzed sequence of {len(sequence)} bp</p>
        <ul>
            <li>Quality Score: <strong>{quality_scores.get('overall', 0)}%</strong></li>
            <li>Detected Repeats: <strong>{len(repeats)}</strong></li>
            <li>High Secondary Structure Regions: <strong>{len(structures)}</strong></li>
        </ul>
        """,
        section_type="text"
    )
    
    # Repeats section
    if repeats:
        repeat_html = "<table><tr><th>Sequence</th><th>Positions</th><th>Copies</th><th>Type</th></tr>"
        for rep in repeats:
            positions_str = ", ".join(map(str, rep.positions[:5]))
            if len(rep.positions) > 5:
                positions_str += f", ... ({len(rep.positions)} total)"
            repeat_type = "Tandem" if rep.spacing == 0 else "Dispersed"
            repeat_html += f"""
            <tr>
                <td><code>{rep.sequence}</code></td>
                <td>{positions_str}</td>
                <td>{rep.copy_count}x</td>
                <td>{repeat_type}</td>
            </tr>
            """
        repeat_html += "</table>"
        report.add_section("Detected Repeats", repeat_html, section_type="table")
    
    # Secondary structures section
    if structures:
        struct_text = "High-stability secondary structure regions detected:\n\n"
        for start, end in structures[:10]:
            struct_text += f"- Position {start}-{end} ({end-start}bp)\n"
        if len(structures) > 10:
            struct_text += f"\n... and {len(structures) - 10} more regions"
        report.add_section("Secondary Structures", struct_text)
    
    return report


def create_optimization_report(original_seq: str, edited_seq: str,
                              metrics: Dict[str, Any]) -> ReportBuilder:
    """
    Create a report comparing original and edited sequences.
    
    Args:
        original_seq: Original sequence
        edited_seq: Optimized sequence
        metrics: Optimization metrics
        
    Returns:
        ReportBuilder with comparison
    """
    from .sequence_utils import gc_percent
    
    report = ReportBuilder("Sequence Optimization Report")
    
    # Key metrics
    report.add_metric("Original Length", len(original_seq), "bp")
    report.add_metric("Edited Length", len(edited_seq), "bp")
    report.add_metric("Original GC%", f"{gc_percent(original_seq):.1f}", "%")
    report.add_metric("Edited GC%", f"{gc_percent(edited_seq):.1f}", "%")
    report.add_metric("Changes Made", metrics.get("num_edits", 0), "edits")
    
    # Summary
    report.add_section(
        "Optimization Summary",
        f"""
        <p>Original sequence: {len(original_seq)} bp</p>
        <p>Edited sequence: {len(edited_seq)} bp</p>
        <p>Number of edits: {metrics.get('num_edits', 0)}</p>
        <p>Motifs removed: {metrics.get('motifs_removed', 0)}</p>
        <p>New motifs introduced: {metrics.get('new_motifs', 0)}</p>
        """,
        section_type="text"
    )
    
    # Edit details
    if "edits" in metrics:
        edit_table = "<table><tr><th>Position</th><th>Original</th><th>Edited</th><th>Type</th></tr>"
        for edit in metrics["edits"][:20]:
            edit_table += f"""
            <tr>
                <td>{edit.get('position', 'N/A')}</td>
                <td><code>{edit.get('original', '')}</code></td>
                <td><code>{edit.get('edited', '')}</code></td>
                <td>{edit.get('type', 'substitution')}</td>
            </tr>
            """
        if len(metrics.get("edits", [])) > 20:
            edit_table += f"<tr><td colspan='4'>... and {len(metrics['edits']) - 20} more edits</td></tr>"
        edit_table += "</table>"
        report.add_section("Edit Details", edit_table, section_type="table")
    
    return report
