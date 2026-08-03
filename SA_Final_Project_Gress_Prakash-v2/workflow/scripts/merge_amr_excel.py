#!/usr/bin/env python3
import sys
import re
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

def merge_to_excel(input_files, output_xlsx):
    # Initialize an Excel writer object
    with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
        # Loop over every input report path
        for path in input_files:
            # Extract sample name from path (e.g. 'sample5') using regex
            match = re.search(r"card_amr_report/([^/]+)/", path)
            sample_name = match.group(1) if match else "Unknown"
            
            # Check if report is empty or exists before reading
            try:
                df = pd.read_csv(path, sep="\t")
            except Exception:
                # Fallback to an empty DataFrame if file fails to load or is empty
                df = pd.DataFrame(columns=["No AMR Genes Detected"])
            
            # Write dataframe to the sheet corresponding to sample_name
            df.to_excel(writer, sheet_name=sample_name, index=False)
            
    # --- Apply Professional Layout and Formatting ---
    wb = openpyxl.load_workbook(output_xlsx)
    
    # Style variables
    font_header = Font(name="Segoe UI", size=11, bold=True, color="FFFFFF")
    font_body = Font(name="Segoe UI", size=10)
    fill_header = PatternFill(start_color="1A365D", end_color="1A365D", fill_type="solid") # Dark Navy
    fill_zebra = PatternFill(start_color="F7FAFC", end_color="F7FAFC", fill_type="solid")  # Light gray-blue
    fill_perfect = PatternFill(start_color="DEF7EC", end_color="DEF7EC", fill_type="solid")# Light green
    fill_strict = PatternFill(start_color="FEF08A", end_color="FEF08A", fill_type="solid") # Light yellow
    
    border_thin = Border(
        left=Side(style='thin', color='E2E8F0'),
        right=Side(style='thin', color='E2E8F0'),
        top=Side(style='thin', color='E2E8F0'),
        bottom=Side(style='thin', color='E2E8F0')
    )
    
    align_center = Alignment(horizontal="center", vertical="center")
    align_left = Alignment(horizontal="left", vertical="center")

    # Iterate back through each worksheet to dress it up
    for ws in wb.worksheets:
        # Explicitly force grid lines to display
        ws.views.sheetView[0].showGridLines = True
        
        # Apply styling column by column, row by row
        for r_idx, row in enumerate(ws.iter_rows(values_only=False), 1):
            for c_idx, cell in enumerate(row, 1):
                cell.font = font_body
                cell.border = border_thin
                
                if r_idx == 1:
                    # Format Header Row
                    cell.font = font_header
                    cell.fill = fill_header
                    cell.alignment = align_center
                else:
                    # Zebra striping for readability
                    if r_idx % 2 == 0:
                        cell.fill = fill_zebra
                    
                    # Data alignments
                    if c_idx in [1, 2, 5, 6]:
                        cell.alignment = align_center
                    else:
                        cell.alignment = align_left
                        
                    # Standardize color-coding for RGI Cut_off classifications
                    col_header = str(ws.cell(row=1, column=c_idx).value or "")
                    if "Cut_Off" in col_header:
                        if cell.value == "Perfect":
                            cell.fill = fill_perfect
                        elif cell.value == "Strict":
                            cell.fill = fill_strict
                            
        # Auto-calculate and adjust column widths
        for col in ws.columns:
            max_len = max(len(str(cell.value or '')) for cell in col)
            col_letter = get_column_letter(col[0].column)
            ws.column_dimensions[col_letter].width = max(max_len + 3, 12)
            
        # Freeze the top header row so it stays pinned on scroll
        ws.freeze_panes = "A2"
        
    wb.save(output_xlsx)

if __name__ == "__main__":
    # Snakemake passes rule arguments via sys.argv:
    # sys.argv[1] is the output file path
    # sys.argv[2:] are all the input file paths
    output_path = snakemake.output.xlsx
    input_paths = snakemake.input.reports
    merge_to_excel(input_paths, output_path)