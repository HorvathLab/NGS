from sklearn.model_selection import train_test_split
from imblearn.over_sampling import RandomOverSampler
from sklearn.preprocessing import StandardScaler
import pandas as pd
from scipy import stats

def pearson(X_train_, X_val_, X_test, y_train_):    
    f = X_train_.apply(lambda c: stats.pearsonr(c,y_train_),axis=0)
    f=f.transpose()
    f.columns = ['rho','pvalue']
    f = pd.DataFrame(f).sort_values(by='rho',ascending=False)
    genes = f.index[f.rho>=.3]
    X_train_ = X_train_.loc[:, X_train_.columns.isin(genes)]
    X_val_ = X_val_.loc[:, X_val_.columns.isin(genes)]
    X_test = X_test.loc[:, X_test.columns.isin(genes)]

    return X_train_, y_train_, X_val_, X_test

def cdk4_6_genes(data_dir, df, genes_file):
    genes = pd.read_csv(data_dir+genes_file, header=None).iloc[:,0].to_list()
    df = df.loc[:, df.columns.isin(genes)]
    
    return df, genes

def cdk4_6_cancer_genes(data_dir, df, cdk4_6_genes_filename, cancer_genes_filename):
    cancer_genes = pd.read_csv(data_dir+cancer_genes_filename,sep='\t')
    cancer_genes = cancer_genes['Gene Symbol'].unique()
    df_, genes = cdk4_6_genes(data_dir, df, cdk4_6_genes_filename)
    df_genes = pd.DataFrame(genes+cancer_genes.tolist())
    genes = df_genes.iloc[:,0].unique().tolist()
    df = df.loc[:, df.columns.isin(genes)]
    
    return df

def split_scale_data(data_dir, df, y_labels, feature_selection, ros=True, cdk4_6_genes_filename=None, cancer_genes_filename=None):
    if feature_selection == 'cdk4_6_genes':
        df, genes = cdk4_6_genes(data_dir, df, cdk4_6_genes_filename)
    elif feature_selection == 'cdk4_6_cancer_genes':
        df = cdk4_6_cancer_genes(data_dir=data_dir, df=df, cdk4_6_genes_filename=cdk4_6_genes_filename, \
            cancer_genes_filename=cancer_genes_filename)

    X_train_, X_valtest, y_train_, y_valtest = train_test_split(df,y_labels,test_size=.4)
    X_val_, X_test, y_val_, y_test = train_test_split(X_valtest,y_valtest, test_size=.5)

    if feature_selection == 'pearson':
        X_train_, y_train_, X_val_, X_test = pearson(X_train_, X_val_, X_test, y_train_)
    
    scaler = StandardScaler()
    X_train_ = scaler.fit_transform(X_train_)
    X_val_ = scaler.transform(X_val_)
    X_test = scaler.transform(X_test)
    if ros:
        ros = RandomOverSampler(random_state=0)
        X_train_ros, y_train_ros = ros.fit_resample(X_train_,y_train_)

    return X_train_ros, y_train_ros, X_val_, y_val_, X_test, y_test