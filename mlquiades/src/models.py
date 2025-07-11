import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn import metrics, tree
from sklearn.metrics import ConfusionMatrixDisplay, auc, f1_score, roc_curve, roc_auc_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeClassifier
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import keras
import keras_tuner
from keras import layers
from plotting import *

def evaluate(y_test, y_pred):
    acc = metrics.accuracy_score(y_test, y_pred)
    f1 = f1_score(np.array(y_test), (y_pred).astype(int))    
    fpr,tpr,thresholds = roc_curve(y_test,y_pred)
    rocauc = roc_auc_score(y_test,y_pred)
    
    return acc, f1.item(), rocauc.item(), fpr, tpr

def evaluate_keras(model, X_test, y_test, y_pred):
    loss, acc = model.evaluate(X_test, y_test)
    fpr, tpr, thresholds_keras = roc_curve(y_test, y_pred)
    rocauc = auc(fpr, tpr)
    f1_metric = keras.metrics.F1Score(threshold=0.5)
    f1_metric.update_state(np.array([y_test.astype(np.float32)]).T,y_pred)
    f1 = f1_metric.result()
    
    return acc, f1_metric, rocauc, fpr, tpr

def decision_tree(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection, plt_confusion=True, plt_rocauc=True):
    clf = tree.DecisionTreeClassifier(max_depth=2)
    clf = clf.fit(X_train_ros, y_train_ros)
    y_pred = clf.predict(X_test)
    
    acc, f1, rocauc, fpr, tpr = evaluate(y_test, y_pred)
    if plt_confusion:
        plot_confusion_matrix(y_test, y_pred, output_dir, feature_selection, model_name='dt', nn=False)
    if plt_rocauc:
        plot_rocauc(fpr, tpr, output_dir, feature_selection, model_name='dt')
    
    return acc, f1, rocauc

def gradient_boosted_decision_tree(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection, plt_confusion=True, plt_rocauc=True):
    clf = GradientBoostingClassifier(n_estimators=2, learning_rate=1.5, max_depth=2,random_state=0).fit(X_train_ros, y_train_ros)
    clf = clf.fit(X_train_ros, y_train_ros)
    y_pred = clf.predict(X_test)
    
    acc, f1, rocauc, fpr, tpr = evaluate(y_test, y_pred)
    if plt_confusion:
        plot_confusion_matrix(y_test, y_pred, output_dir, feature_selection, model_name='gbdt', nn=False)
    if plt_rocauc:
        plot_rocauc(fpr, tpr, output_dir, feature_selection, model_name='gbdt')
    
    return acc, f1, rocauc

def neural_net(X_train_ros, y_train_ros, X_val_, y_val_, X_test, y_test, output_dir, feature_selection, plt_confusion=True, plt_rocauc=True):
    model = Sequential()
    model.add(Dense(12, input_shape=(X_train_ros.shape[1],), activation='relu'))
    model.add(Dense(8, activation='relu'))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    
    callback = keras.callbacks.EarlyStopping(monitor='val_loss',patience=3,min_delta=.01)
    optimiizer = keras.optimizers.Adam(learning_rate=.0001)

    model.compile(loss='binary_crossentropy', optimizer=optimiizer, metrics=['accuracy'])
    y_train_ros = y_train_ros.replace(-1,0)
    y_val_ = y_val_.replace(-1,0)
    model.fit(X_train_ros, y_train_ros, validation_data = (X_val_, y_val_), epochs=100, batch_size=50, callbacks=[callback])
    y_pred = model.predict(X_test)
    y_test = y_test.replace(-1,0)
    
    acc, f1, rocauc, fpr, tpr = evaluate_keras(model, X_test, y_test, y_pred)
    if plt_confusion:
        plot_confusion_matrix(y_test, y_pred, output_dir, feature_selection, model_name='nn', nn=True)
    if plt_rocauc:
        plot_rocauc(fpr, tpr, output_dir, feature_selection, model_name='nn')
    
    return acc, f1, rocauc

def neural_net_with_hyperband(X_train_ros, y_train_ros, X_val_, y_val_, X_test, y_test, data_dir, step_size_nodes, 
                              min_nodes, max_nodes, max_trials, executions_per_trial, patience, min_delta,
                              epochs, learning_rate_min, learning_rate_max, output_dir, feature_selection, 
                              plt_confusion=True, plt_rocauc=True):
    def build_model(hp):
        model = keras.Sequential()
        for i in range(hp.Int('n_layers', min_value=1, max_value=20)):
            model.add(Dense(units=hp.Int('units__' + str(i),
                                        min_value=min_nodes,
                                        max_value=max_nodes,
                                        step=step_size_nodes),
                            activation='relu'))
            # model.add(Dropout(hp.Float('drpout', min_value=0, max_value=.9, step=.1)))
        model.add(layers.Dense(1, activation='sigmoid'))
        learning_rate = hp.Float('learning_rt', min_value=learning_rate_min, max_value=learning_rate_max, sampling='log')
        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
        return model

    build_model(keras_tuner.HyperParameters())

    tuner = keras_tuner.RandomSearch(
        hypermodel=build_model,
        objective='val_accuracy',
        max_trials=max_trials,
        executions_per_trial=executions_per_trial,
        overwrite=True,
        directory=data_dir
    )

    tuner.search_space_summary()
    
    stop_early = keras.callbacks.EarlyStopping(monitor='val_loss', patience=patience, min_delta=min_delta)
    tuner.search(X_train_ros, y_train_ros.replace(-1,0), epochs=epochs, validation_data=(X_val_, y_val_), callbacks=[stop_early])
    best_model = tuner.get_best_models(num_models=1)[0]
    y_pred = best_model.predict(X_test)
    
    acc, f1, rocauc, fpr, tpr = evaluate_keras(best_model, X_test, y_test, y_pred)
    if plt_confusion:
        plot_confusion_matrix(y_test, y_pred, output_dir, feature_selection, model_name='nn_hb', nn=True)
    if plt_rocauc:
        plot_rocauc(fpr, tpr, output_dir, feature_selection, model_name='nn_hb')
    
    return acc, f1, rocauc

def random_forest(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection, max_depth=2, random_state=0, \
    plt_confusion=True, plt_rocauc=True):
    clf = RandomForestClassifier(max_depth=max_depth, random_state=random_state).fit(X_train_ros, y_train_ros)
    y_pred = clf.predict(X_test)
    
    acc, f1, rocauc, fpr, tpr = evaluate(y_test, y_pred)
    if plt_confusion:
        plot_confusion_matrix(y_test, y_pred, output_dir, feature_selection, model_name='rf', nn=False)
    if plt_rocauc:
        plot_rocauc(fpr, tpr, output_dir, feature_selection, model_name='rf')
    
    return acc, f1, rocauc

def ridge_classifier(X_train_ros, y_train_ros, X_test, y_test, output_dir, feature_selection, plt_confusion=True, plt_rocauc=True):
    clf = RidgeClassifier().fit(X_train_ros, y_train_ros)
    y_pred = clf.predict(X_test)
    
    acc, f1, rocauc, fpr, tpr = evaluate(y_test, y_pred)
    if plt_confusion:
        plot_confusion_matrix(y_test, y_pred, output_dir, feature_selection, model_name='ridge', nn=False)
    if plt_rocauc:
        plot_rocauc(fpr, tpr, output_dir, feature_selection, model_name='ridge')
    
    return acc, f1, rocauc
