#!/usr/bin/env python3
"""
Convert rsIDs to hs37d5 genomic coordinates
hs37d5 is the 1000 Genomes Project reference (GRCh37 + decoy sequences)
"""

import requests
import json
import time
import sys


def get_hs37d5_loci_myvariant(rsid):
    """
    Get hs37d5 coordinates using MyVariant.info API.
    hs37d5 uses the same coordinates as GRCh37/hg19 for chromosomes 1-22, X, Y, MT.
    """
    if not rsid.startswith('rs'):
        rsid = f"rs{rsid}"
    
    url = f"https://myvariant.info/v1/variant/{rsid}"
    params = {
        'assembly': 'hg19',  # hg19 and GRCh37 have same coordinates as hs37d5 main chromosomes
        'fields': 'dbsnp.rsid,dbsnp.chrom,dbsnp.hg19,dbsnp.ref,dbsnp.alt'
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        
        if response.status_code == 404:
            return None
        
        response.raise_for_status()
        data = response.json()
        
        result = {
            'rsid': rsid,
            'chromosome': None,
            'position': None,
            'ref_allele': None,
            'alt_alleles': [],
            'build': 'hs37d5'
        }
        
        # Extract coordinates
        if 'dbsnp' in data and 'hg19' in data['dbsnp']:
            hg19_data = data['dbsnp']['hg19']
            
            # Handle both dict and list formats
            if isinstance(hg19_data, dict):
                chrom = str(hg19_data.get('chr', '')).replace('chr', '')
                result['chromosome'] = chrom
                result['position'] = hg19_data.get('start')
            elif isinstance(hg19_data, list) and len(hg19_data) > 0:
                chrom = str(hg19_data[0].get('chr', '')).replace('chr', '')
                result['chromosome'] = chrom
                result['position'] = hg19_data[0].get('start')
            
            # Get alleles
            if 'ref' in data['dbsnp']:
                result['ref_allele'] = data['dbsnp']['ref']
            if 'alt' in data['dbsnp']:
                alt = data['dbsnp']['alt']
                result['alt_alleles'] = [alt] if isinstance(alt, str) else alt
        
        if result['chromosome'] and result['position']:
            return result
        return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error with MyVariant.info API: {e}")
        return None


def get_hs37d5_loci_ensembl(rsid):
    """
    Get hs37d5 coordinates using Ensembl REST API (GRCh37 archive).
    hs37d5 main chromosomes match GRCh37 coordinates.
    """
    if not rsid.startswith('rs'):
        rsid = f"rs{rsid}"
    
    url = f"https://grch37.rest.ensembl.org/variation/human/{rsid}"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(url, headers=headers, timeout=10)
        
        if response.status_code == 404:
            return None
        
        response.raise_for_status()
        data = response.json()
        
        result = {
            'rsid': rsid,
            'chromosome': None,
            'position': None,
            'ref_allele': None,
            'alt_alleles': [],
            'build': 'hs37d5'
        }
        
        if 'mappings' in data and len(data['mappings']) > 0:
            mapping = data['mappings'][0]
            
            chrom = mapping.get('seq_region_name', '').replace('chr', '')
            result['chromosome'] = chrom
            result['position'] = mapping.get('start')
            
            if 'allele_string' in mapping:
                alleles = mapping['allele_string'].split('/')
                if len(alleles) > 0:
                    result['ref_allele'] = alleles[0]
                if len(alleles) > 1:
                    result['alt_alleles'] = alleles[1:]
        
        if result['chromosome'] and result['position']:
            return result
        return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error with Ensembl API: {e}")
        return None


def get_hs37d5_loci(rsid, method='myvariant'):
    """
    Get hs37d5 genomic coordinates for a given rsID.
    
    Args:
        rsid: SNP rsID (with or without 'rs' prefix)
        method: 'myvariant', 'ensembl', or 'auto' for fallback
    
    Returns:
        dict with chromosome, position, and other info, or None if not found
    
    Note:
        hs37d5 is the 1000 Genomes Project reference genome, which is:
        - GRCh37 primary assembly (chromosomes 1-22, X, Y, MT)
        - Plus decoy sequences (hs37d5 specific)
        - Plus unlocalized/unplaced contigs
        
        For standard chromosomes, coordinates match GRCh37/hg19 exactly.
    """
    if not rsid.startswith('rs'):
        rsid = f"rs{rsid}"
    
    if method == 'myvariant':
        return get_hs37d5_loci_myvariant(rsid)
    elif method == 'ensembl':
        return get_hs37d5_loci_ensembl(rsid)
    elif method == 'auto':
        # Try MyVariant first, fall back to Ensembl
        result = get_hs37d5_loci_myvariant(rsid)
        if result is None:
            result = get_hs37d5_loci_ensembl(rsid)
        return result
    else:
        raise ValueError(f"Unknown method: {method}")


def batch_get_loci(rsid_list, delay=0.2, method='myvariant'):
    """
    Get hs37d5 coordinates for multiple rsIDs.
    
    Args:
        rsid_list: List of rsIDs
        delay: Delay between requests in seconds
        method: API method to use
    
    Returns:
        dict mapping rsIDs to their loci information
    """
    results = {}
    
    for i, rsid in enumerate(rsid_list):
        if not rsid.startswith('rs'):
            rsid = f"rs{rsid}"
            
        result = get_hs37d5_loci(rsid, method=method)
        if result:
            results[rsid] = result
        else:
            print(f"Could not find coordinates for {rsid}")
        
        # Rate limiting
        if delay > 0 and i < len(rsid_list) - 1:
            time.sleep(delay)
    
    return results


def format_loci(loci, style='standard'):
    """
    Format loci information for display.
    
    Args:
        loci: Loci dictionary
        style: 'standard' (chr1:12345) or 'hs37d5' (1:12345 - no 'chr' prefix)
    """
    if not loci:
        return "Not found"
    
    chrom = loci['chromosome']
    pos = loci['position']
    ref = loci.get('ref_allele', 'N/A')
    alt = ','.join(loci.get('alt_alleles', [])) or 'N/A'
    
    if style == 'hs37d5':
        # hs37d5 typically uses no 'chr' prefix
        return f"{chrom}:{pos} ({ref}/{alt})"
    else:
        # Standard format with 'chr' prefix
        return f"chr{chrom}:{pos} ({ref}/{alt})"


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Command line usage
        rsids = sys.argv[1:]
        print(f"Looking up {len(rsids)} rsID(s) for hs37d5...\n")
        
        for rsid in rsids:
            loci = get_hs37d5_loci(rsid, method='auto')
            if loci:
                # Use hs37d5 style (no chr prefix) for display
                print(f"{loci['rsid']}: {format_loci(loci, style='hs37d5')}")
            else:
                print(f"{rsid}: Not found")
    else:
        # Example usage
        print("hs37d5 Assembly Coordinate Lookup")
        print("=" * 60)
        print("Note: hs37d5 is the 1000 Genomes Project reference")
        print("      Coordinates for chr 1-22,X,Y,MT match GRCh37/hg19")
        print("=" * 60)
        print("\nSingle rsID lookup (MyVariant.info):")
        print("-" * 60)
        
        rsid = "rs429358"  # APOE SNP
        loci = get_hs37d5_loci(rsid, method='myvariant')
        
        if loci:
            print(f"rsID: {loci['rsid']}")
            print(f"Chromosome: {loci['chromosome']}")
            print(f"Position: {loci['position']}")
            print(f"Locus (hs37d5 style): {format_loci(loci, style='hs37d5')}")
            print(f"Locus (standard): {format_loci(loci, style='standard')}")
            print(f"Ref allele: {loci['ref_allele']}")
            print(f"Alt alleles: {', '.join(loci['alt_alleles'])}")
        
        print("\n" + "="*60)
        print("Batch lookup:")
        print("-" * 60)
        
        # Multiple rsIDs
        rsids = ["rs7412", "rs429358", "rs1799945", "rs1800562"]
        results = batch_get_loci(rsids, delay=0.2, method='auto')
        
        for rsid, loci in results.items():
            if loci:
                print(f"{loci['rsid']}: {format_loci(loci, style='hs37d5')}")
        
        print("\n" + "="*60)
        print("\nCommand line usage:")
        print("  python3 rsid_to_hs37d5.py rs429358 rs7412 rs1799945")
        print("\nFor IGV with hs37d5 reference:")
        print("  - Use chromosome names without 'chr' prefix (e.g., '1' not 'chr1')")
        print("  - Coordinates are identical to GRCh37/hg19 for standard chromosomes")
