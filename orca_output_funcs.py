#!/usr/bin/env python
# coding: utf-8

# In[229]:


import os
import csv

    
def extract_data_from_orca_output(directory = '/home/biswajit/Documents/Naman/analysis/sp', output_csv_name = 'output_csv', \
                            file_name='CN_6', file_type ='SP'):
    data_list = []
    d = []
    
    def is_float(value):
        try:
            float(value)
            return True
        except ValueError:
            return False
    
    
    if file_type == 'Opt':
        # Geo opt out file
        search_keywords = ['Final Gibbs free energy', 'G-E(el)']
        freq_keyword = ['VIBRATIONAL FREQUENCIES']
        f_name  = directory + '/' + '%s_opt.out'%(file_name)
    elif file_type == 'SP':
        # Single point out file
        search_keywords = ['FINAL SINGLE POINT ENERGY']
        f_name  = directory + '/' + '%s_sp.out'%(file_name)
    
    
    
    with open(f_name, 'r') as infile:
        # check for negative frequency
        freq_line_ndx = []
        if file_type == 'Opt':
            print('sss')
            for line_number, line in enumerate(infile, start=1):
                #if any(keyword in line for keyword in freq_keyword):
                if line.startswith(freq_keyword[0]):

                    freq_line_ndx.append(line_number)
                elif line.startswith('NORMAL MODES'):
                    freq_line_ndx.append(line_number)
                    
        
        infile.seek(0)
        lines = infile.readlines()
        freq_lines = lines[freq_line_ndx[0]+4:freq_line_ndx[1]-4]
        fl = []
        for v in freq_lines:
            v = [float(word) for word in v.split() if is_float(word)]

            fl.append(v[0])
        negative_frequencies = [v for v in fl if v < 0]
        if len(negative_frequencies) > 1:
            print("WARNING: imaginary frequencies ({}) detected in geometry optimization file of {}".format(negative_frequencies, f_name))
            
        else:
            pass
        
        # get energy
        infile.seek(0)
        for line_number, line in enumerate(infile, start=1):
            if any(keyword in line for keyword in search_keywords):
                # Extract float values from the line
                float_values = [float(word) for word in line.split() if is_float(word)]
                #data_list[keyword] = float_values
                d.append(float_values)
                # Append data to the list
                data_list.append({
                    'Filename': f_name,
                    'Line': line_number,
                    'FloatValues': float_values
                })

    
    if file_type == 'Opt':
        energy_hartree = d[0][0]
        free_energy_correction = d[1][0]
        return energy_hartree, free_energy_correction, data_list
    
    else:
        energy_hartree = d[0]
        return energy_hartree, data_list


# In[230]:


extract_data_from_orca_output(directory = '/home/biswajit/Documents/Naman/analysis/sp', output_csv_name = 'output_csv', \
                            file_name='CN_6', file_type ='Opt')


# In[184]:


def get_bfe_from_orca(reactant = ['UO2_5water', 'P2A'], reactant_unit = [1,2], \
            product = ['CN_5', 'water_4'], product_unit = [1, 1], \
            directory = '/home/biswajit/Documents/Naman/analysis/sp', corrected_bfe = True, unit = 'kJ/mol'):
    
    
    """
    reactant: file name of orca-generated output files for the reactants. mention name without _sp.out or _opt.out
    product: file name of orca-generated output files for the products. mention name without _sp.out or _opt.out
    
    NOTE: 
    the input file name must end with '_sp' and '_opt' for single point and geometry optimization 
    output files, respectively.
    
    reactant_unit = the stoichiometry (of reactants) for the reaction in array form.
    product_unit = the stoichiometry (of products) for the reaction in array form.
    
    corrected_bfe = if set to True, free energy corrected binding energy will be computed. 
    Please ensure geometry opt files are present in the same directory for this computations
    
    unit: unit of binding free energy. either kcal/mol or kJ/mol 
    
    
    """
    energies = []
    
    if corrected_bfe:
        for each in (reactant + product):
            energy_hartree, _ = extract_data_from_orca_output(directory = directory, \
                                file_name=each, file_type ='SP')

            e, correction, d = extract_data_from_orca_output(directory = directory, \
                                file_name=each, file_type ='Opt')
            
            energy_hartree = energy_hartree[0] + correction


            energies.append([energy_hartree])
    else:
        
        for each in (reactant + product):
            energy_hartree, _ = extract_data_from_orca_output(directory = '/home/biswajit/Documents/Naman/analysis/sp', output_csv_name = 'output_csv', \
                                file_name=each, file_type ='SP')

            energies.append(energy_hartree)
        
    all_sp = reactant + product
    
    r_sum = 0
    for _ in range(len(reactant_unit)):
        r_sum += reactant_unit[_]*energies[_][0]
    
    p_sum = 0
    for _ in range(len(product_unit)):
        val = len(reactant_unit) + _
        p_sum += product_unit[_] * energies[val][0]
        
    if unit == 'kJ/mol':  
        bfe = ((p_sum - r_sum)*627.5095)*4.18
    elif unit == 'kcal/mol':
        bfe = ((p_sum - r_sum)*627.5095)
          
    return bfe


# In[185]:


get_bfe_from_orca(reactant = ['UO2_5water', 'P2A'], reactant_unit = [1,2], \
            product = ['CN_4', 'water_5'], product_unit = [1,1], \
            directory = '/home/biswajit/Documents/Naman/analysis/sp', corrected_bfe = True, unit = 'kcal/mol')


# In[ ]:




