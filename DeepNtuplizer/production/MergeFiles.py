import uproot
import pandas as pd
import os
from glob import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

def process_file(file_path, tree_name, chunk_size=100000):
    """Process a single ROOT file in chunks with event selection"""
    selected_chunks = []
    
    try:
        # Open the file and get the tree
        with uproot.open(file_path) as file:
            if tree_name not in file:
                print(f"Tree {tree_name} not found in {file_path}")
                return pd.DataFrame()

            # Process the tree in chunks
            for chunk in uproot.iterate(
                file[tree_name], 
                step_size=chunk_size, 
                library="pd"
            ):
                # Apply selection criteria
                filtered = select_events(chunk)
                if not filtered.empty:
                    selected_chunks.append(filtered)

    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")
        return pd.DataFrame()

    return pd.concat(selected_chunks, ignore_index=True) if selected_chunks else pd.DataFrame()

def process_root_files_parallel(
    folders_sgn,
    folders_bkg, 
    output_file="combined_events.root",
    tree_name="myTree",
    max_workers=None,
    chunk_size=100000
):
    """
    Process ROOT files from multiple folders in parallel with chunked reading
    
    Args:
        folders (list): List of paths to folders containing ROOT files
        output_file (str): Output ROOT file name
        tree_name (str): Name of TTree in ROOT files
        max_workers (int): Maximum number of parallel workers
        chunk_size (int): Number of entries per chunk
    """
    
    # Collect all ROOT files from all folders
    all_files_sgn = []
    for folder in folders_sgn:
        file_pattern = os.path.join(folder, "*.root")
        root_files = glob(file_pattern)
        all_files_sgn.extend(root_files)
   # print("sgn", all_files_sgn)
    
    all_files_bkg = []
    for folder in folders_bkg:
        file_pattern = os.path.join(folder, "*.root")
        root_files = glob(file_pattern)
        all_files_bkg.extend(root_files)

   # print("bkg", all_files_bkg)
    
    if not all_files_sgn:
        print("No ROOT files found in specified folders")
        return pd.DataFrame()
    if not all_files_bkg:
        print("No ROOT files found in specified folders")
        return pd.DataFrame()

    
    # Create partial function for fixed parameters
    process_file_partial = partial(
        process_file,
        tree_name=tree_name,
        chunk_size=chunk_size
    )

    combined_df = pd.DataFrame()
    '''
    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_file_partial, file): file for file in all_files}
        
        for future in as_completed(futures):
            file_path = futures[future]
            try:
                result = future.result()
                if not result.empty:
                    combined_df = pd.concat([combined_df, result], ignore_index=True)
                    print(f"Processed {len(result)} events from {file_path}")
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")

    # Save combined events to new ROOT file
    if not combined_df.empty:
        with uproot.recreate(output_file) as f:
            f[tree_name] = combined_df
        print(f"Combined {len(combined_df)} events saved to {output_file}")
    else:
        print("No events selected from any files")
    
    return combined_df
    '''
def select_signal(df):
    """Example selection criteria - modify for your analysis"""
    # Example: Select events with energy > 0.5 GeV and status == 1
    mask = df['isB'] or df['isBB'] or df['isGBB'] or df['isLeptoncB'] or df['isLeptonicB_C'] or df['isC'] or df['isCC'] or df['isGCC']
    return df[mask]

def select_bkg(df):
    """Example selection criteria - modify for your analysis"""
    # Example: Select events with energy > 0.5 GeV and status == 1
    mask = df['isU'] or df['isD'] or df['isS'] or df['isG']
    return df[mask]

if __name__ == "__main__":
    # List of folders containing ROOT files
    folders_bkg = [
         "/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v2/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/PFC_QCD_pt15to7k_ext1/250326_153555/0000/"
    ]
    folders_sgn = [
         "/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v3/CRAB_UserFiles/PFC_Signal_chain_m70_dm20_13_04_25/250413_155615/0000/"
    ]
    
    # Create output directory if needed
    os.makedirs("output", exist_ok=True)
    
    # Process files with parallel execution
    combined_data = process_root_files_parallel(
        folders_sgn,
        folders_bkg,
        output_file="output/combined_events.root",
        tree_name="deepntuplizer/tree", 
        max_workers=os.cpu_count(),  
        chunk_size=100000  # memory constraints
    )
