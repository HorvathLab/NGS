import os
import argparse
from preprocessing import *
from processing import *
from models import *
from plotting import *

def params():
    parser = argparse.ArgumentParser()
    parser.add_argument('--a', type=str, action='store', dest='data_dir', 
                        help='data directory')
    parser.add_argument('--b', type=str, action='store', dest='output_folder_name', default='output',
                        help='output folder name')
    parser.add_argument('--c', type=str, action='store', dest='ccle_file', 
                        help='ccle data file')
    parser.add_argument('--d', type=str, action='store', dest='drug_file', 
                        help='drug data file')
    parser.add_argument('--e', type=int, action='store', dest='ic50_cutoff_value', 
                        help='ic50 cutoff value')
    parser.add_argument('--f', type=bool, action='store', dest='ros', default=True,
                        help='randomly oversample training data')
    parser.add_argument('--g', type=int, action='store', dest='step_size_nodes', default=5,
                        help='step size for nodes in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--i', type=int, action='store', dest='min_nodes', default=5,
                        help='minimum number of nodes in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--j', type=int, action='store', dest='max_nodes', default=512,
                        help='maximum number of nodes in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--k', type=int, action='store', dest='max_trials', default=10,
                        help='maximum number of trials in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--l', type=int, action='store', dest='executions_per_trial', default=3,
                        help='number of executions per trial in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--m', type=int, action='store', dest='patience', default=3,
                        help='patience in early stopping for nn w/ hyperband')
    parser.add_argument('--n', type=int, action='store', dest='min_delta', default=.01,
                        help='minimum delta in early stopping for nn w/ hyperband')
    parser.add_argument('--o', type=int, action='store', dest='epochs', default=100,
                        help='epochs in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--p', type=int, action='store', dest='learning_rate_min', default=1e-4,
                        help='learning rate minimum in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--q', type=int, action='store', dest='learning_rate_max', default=1e-2,
                        help='learning rate maximum in hyperparameter tuning for nn w/ hyperband')
    parser.add_argument('--r', type=str, action='store', dest='feature_selection',
                        default=None, help='choose one of three options for feature selection \
                        for the ml models, options include: cdk4_6_genes, cdk4_6_cancer_genes, pearson')
    parser.add_argument('--s', type=str, action='store', dest='cdk4_6_genes_filename',
                        help='filename for cdk4 and cdk6 genes for the feature selection option cdk4_6_genes')
    parser.add_argument('--t', type=str, action='store', dest='cancer_genes_filename',
                        help='filename for cancer genes for the feature selection option cdk4_6_cancer_genes')
    parser.add_argument('--u', type=str, action='store', dest='genes_gtf',
                        help='filename for genes.gtf file (e.g. gencode.v19.genes.v7_model.patched_contigs.gtf)')
    return parser

def main():
    parser = params()
    args = parser.parse_args()
    data_dir = args.data_dir+'/'
    output_dir = args.output_folder_name
    ccle_file = args.ccle_file
    drug_file = args.drug_file
    ic50_cutoff_value = args.ic50_cutoff_value
    ros = args.ros
    step_size_nodes = args.step_size_nodes
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    max_trials = args.max_trials
    executions_per_trial = args.executions_per_trial
    patience = args.patience
    min_delta = args.min_delta
    epochs = args.epochs
    learning_rate_min = args.learning_rate_min
    learning_rate_max = args.learning_rate_max
    feature_selection = args.feature_selection
    cdk4_6_filename = args.cdk4_6_genes_filename
    cancer_genes_filename = args.cancer_genes_filename
    genes_gtf = args.genes_gtf
    
    if feature_selection is None:
        raise TypeError('missing feature selection (option --r)')
    if genes_gtf is None:
        raise TypeError('missing genes.gtf (option --u)')
    if feature_selection == 'cdk4_6_genes':
        if cdk4_6_filename is None:
            raise TypeError('missing cdk4_6_genes_filename (option --s)')
    if feature_selection == 'cdk4_6_cancer_genes':
        if cancer_genes_filename is None:
            raise TypeError('missing cancer_genes_filename (option --t)')

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    df, y_labels = read_in_data(data_dir, ccle_file, drug_file, ic50_cutoff_value, genes_gtf)
    X_train_ros, y_train_ros, X_val_, y_val_, X_test, y_test = split_scale_data(data_dir=data_dir, \
        df=df, y_labels=y_labels, ros=ros, feature_selection=feature_selection, \
        cdk4_6_genes_filename=cdk4_6_filename, cancer_genes_filename=cancer_genes_filename)
    
    evaluation_df = []
    evaluation_df.append(decision_tree(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection,))
    evaluation_df.append(gradient_boosted_decision_tree(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection,))
    evaluation_df.append(neural_net(X_train_ros, y_train_ros, X_val_, y_val_, X_test, y_test, output_dir, feature_selection,))
    evaluation_df.append(neural_net_with_hyperband(X_train_ros, y_train_ros, X_val_, y_val_, X_test, 
                                                   y_test, data_dir, step_size_nodes, min_nodes, 
                                                   max_nodes, max_trials, executions_per_trial, 
                                                   patience, min_delta, epochs, learning_rate_min, 
                                                   learning_rate_max, output_dir, feature_selection,))
    evaluation_df.append(random_forest(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection))
    evaluation_df.append(ridge_classifier(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection,))
    plot_combined_rocauc(evaluation_df, feature_selection, output_dir)

if __name__=='__main__':
    main()