import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import ConfusionMatrixDisplay

def plot_combined_rocauc(evaluation_df, feature_selection, output_dir):
    evaluation_df = pd.DataFrame(evaluation_df)
    evaluation_df.columns = ['acc','f1', 'rocauc']
    evaluation_df['model'] = ['dt','gbdt','nn','nn w/ hyperband','random forest','ridgeclassifier']

    df__ = evaluation_df.drop(columns=['acc','f1']).melt('model')
    df__['value'] = np.array(df__['value'],dtype='float')
    df__ = pd.DataFrame({'model': df__['model'], 'variable': df__['variable'],'value': df__['value']})

    bp = sns.barplot(df__, x='variable', y='value', hue='model')
    bp.set(xlabel=' ', ylabel=' ', title=feature_selection)
    fig = bp.get_figure()
    fig.savefig(output_dir+'/'+'plt_rocauc_all_'+feature_selection+'.png')

def plot_confusion_matrix(y_test, y_pred, output_dir, feature_selection, model_name, nn=False):
    if not nn:
        ConfusionMatrixDisplay.from_predictions(np.array([y_test.astype(np.float32)]).T,(y_pred).astype(int),display_labels=[0,1])
        plt.savefig(output_dir+'/'+'plt_confusion_'+model_name+'_'+feature_selection+'.png')
        plt.close()
    else:
        ConfusionMatrixDisplay.from_predictions(np.array([y_test.replace(-1,0).astype(np.float32)]).T,(y_pred>=.5).astype(int),display_labels=[0,1])
        plt.savefig(output_dir+'/'+'plt_confusion_'+model_name+'_'+feature_selection+'.png')
        plt.close()

def plot_rocauc(fpr, tpr, output_dir, feature_selection, model_name):
    plt.plot(fpr,tpr)
    plt.plot([0,1],[0,1],'--')
    plt.title('ROC AUC '+model_name+' '+feature_selection)
    plt.xlabel('fpr')
    plt.ylabel('tpr')
    plt.savefig(output_dir+'/'+'plt_rocauc_'+model_name+'_'+feature_selection+'.png')
    plt.close()

