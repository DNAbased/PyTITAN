def main():
    import argparse
    import pandas as pd
    import numpy as np
    import logging
    import concurrent.futures
    from handlers import ConversionHandler, DataHandler, BindingHandler
    
    parser = argparse.ArgumentParser(description='Analyze the impact of genomic variants on transcription factor binding affinity', usage='%(prog)s [options]')
    parser.add_argument('-i', '--input', type=str, metavar='input', help='Input file of variants to be analyzed. Tab-separated. Assumes that a header is present (required columns: CHROM, POS, REF, ALT).')
    parser.add_argument('-f', '--folder', type=str, metavar='folder', help='Folder containing the transcription factor binding site information.')
    parser.add_argument('-o', '--output', type=str, default='out.tsv', metavar='output', help='Output file path.')
    parser.add_argument('-g', '--genome', type=str, default='hs37d5.fa', metavar='genome', help='Path to genome in FASTA format.')
    parser.add_argument('-t', '--threshold', type=int, default=50, metavar='threshold', help='Threshold (percentile) for binding affinity of a variant position / transcription factor binding site combination to be included.')
    parser.add_argument('-m', '--min_icp', type=float, default=1.0, metavar='min_icp', help='Minimum information content of position.') # e.g. Minimum information content of affected positions to be reported (0 to 2). # Thus filter on total information content of position
    parser.add_argument('-n', '--min_icb', type=float, default=0.5, metavar='min_icb', help='Minimum information content of base.')
    parser.add_argument('-p', '--processes', type=int, default=1, metavar='processes', help='Number of processes to use.')
    parser.add_argument('-q', '--pseudocount', type=float, default=0.8, metavar='pseudocount', help='Pseudocount to use for PWM generation.')
    parser.add_argument('-b', '--background', type=int, nargs=4, default=[0.25, 0.25, 0.25, 0.25], metavar='background', help='Nucleotide background probabilities for ACGT.')
    args = parser.parse_args()
    
    # load handlers
    ch = ConversionHandler()
    dh = DataHandler()
    bh = BindingHandler()

    # set up logger
    logger = logging.getLogger(__name__)
    log_handler = logging.FileHandler(args.output + '.log')
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    log_format = logging.Formatter('%(asctime)s.%(msecs)03d - %(levelname)s: %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    log_handler.setFormatter(log_format)
    logger.addHandler(log_handler)

    nucleotides = ['A', 'C', 'G', 'T']
    bg_prob = dict(zip(nucleotides, args.background))

    logger.info('---Analysis started---')
    logger.info('Input file: ' + str(args.input))
    logger.info('Matrix folder: ' + str(args.folder))
    logger.info('Output path: ' + str(args.output))
    logger.info('Genome path: '+ str(args.genome))
    logger.info('Percentile threshold: ' + str(args.threshold))
    logger.info('Minimum position IC: ' + str(args.min_icp))
    logger.info('Minimum base IC: ' + str(args.min_icb))
    logger.info('Number of processes: ' + str(args.processes))
    logger.info('Pseudocount: ' + str(args.pseudocount))
    logger.info('Background probabilities: ' + str(bg_prob))

    folder_input = dh.get_files(filepath=args.folder)
    logger.info('Number of matrices: ' + str(len(folder_input)))
    
    variants = pd.read_csv(args.input, sep='\t')
    logger.info('Number of input variants: ' + str(len(variants)))

    variants = variants.drop_duplicates()
    variants = variants.reset_index(drop=True)

    if 'ID' not in variants.columns:
        variants = dh.create_id(df=variants)
    
    variants = dh.get_seq(df=variants, genome=args.genome)
    variants = dh.alt_seq(df=variants)
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        spawned_processes = executor.map(bh.check_variants, [variants]*args.processes, np.array_split(folder_input, args.processes), [args.threshold]*args.processes, [args.pseudocount]*args.processes, [bg_prob]*args.processes, [args.min_icp]*args.processes, [args.min_icb]*args.processes) # possible to name the arguments here?
    
    variant_result = pd.concat(spawned_processes, sort=False)

    final_result = bh.filter_variants(df=variant_result)
    
    final_result.to_csv(args.output, sep='\t', index=False)
    logger.info('---Analysis finished---')

if __name__ == '__main__':
    main()

# for every tfm, check sequence, put score into list, move one base forward, repeat, finally get max value from list
# for minus strand, the position of the variant is reversed, i.e. position 8/8 is actually 1/8 because it is only 8/8 for the reverse complemented variant sequence
