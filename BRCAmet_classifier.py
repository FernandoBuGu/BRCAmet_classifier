#Project: Deep Learning Classifier to estimate the probability of future metastasis from Breast Primary Tumor WES data
#@feBueno fernando.bueno.gutierrez@usal.es August-2020

#INPUT:
#Six .txt generated with BRCAmet_getInput.R: Xtr.txt, Xv.txt, Xte.txt, ytr.txt, yv.txt, yte.txt
#with predictors and label binary matrixes data for train, validation and test sets repectivelly

#OUTPUT:
#test accuracy (80%) on a balanced set with 100 samples from each class


import tensorflow as tf
from tensorflow import keras
import tensorflow_datasets as tfds
import numpy as np
import json; import math; import random; import os; import csv; import zipfile#; import sklearn
import matplotlib.pyplot as plt

physical_devices = tf.config.experimental.list_physical_devices('GPU')
assert len(physical_devices) > 0, "Not enough GPU hardware devices available"
config = tf.config.experimental.set_memory_growth(physical_devices[0], True)

X_train = np.loadtxt("/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/Xtr.txt",skiprows=0)
X_valid = np.loadtxt("/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/Xv.txt",skiprows=0)
y_train = np.loadtxt("/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/ytr.txt",skiprows=0)
y_valid = np.loadtxt("/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/yv.txt",skiprows=0)

model = tf.keras.models.Sequential([
    tf.keras.layers.Dense(50, activation=tf.nn.relu),
    tf.keras.layers.Dropout(0.8),
    tf.keras.layers.Dense(50, activation=tf.nn.relu),
    tf.keras.layers.Dropout(0.8),
    tf.keras.layers.Dense(1, activation=tf.nn.sigmoid)
])

moptimizer=tf.keras.optimizers.Adam()#keras.optimizers.SGD(lr=0.0001)
model.compile(optimizer=moptimizer,loss='binary_crossentropy',metrics=['accuracy'])

cp_filepath = '/home/febueno/PycharmProjects/BRCAmet/cp.h5'
model_checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
    filepath=cp_filepath,
    save_weights_only=False,
    monitor='val_accuracy',
    mode='auto',
    save_best_only=True,
    verbose=1)

history = model.fit(X_train, y_train, epochs=300, validation_data=(X_valid, y_valid), batch_size=32,callbacks=[model_checkpoint_callback])


plt.figure(figsize=(10, 6))
def plot_graphs(history, string):
  plt.plot(history.history[string])
  plt.plot(history.history['val_'+string])
  plt.xlabel("Epochs")
  plt.ylabel(string)
  plt.legend([string, 'val_'+string])
  plt.show()
plot_graphs(history, "accuracy")
plot_graphs(history, "loss")
plt.legend(["accuracy", 'val_'+"accuracy","loss", 'val_'+"loss"])


x_test = np.loadtxt("/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/Xte.txt",skiprows=0)
y_test = np.loadtxt("/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/yte.txt",skiprows=0)

te_unique, te_counts = np.unique(y_test, return_counts=True)
te_D=dict(zip(te_unique, te_counts))
print(te_D)

#load best model
    #restart model
tf.keras.backend.clear_session()
model = tf.keras.models.Sequential([
    tf.keras.layers.Dense(50, activation=tf.nn.relu),
    tf.keras.layers.Dropout(0.8),
    tf.keras.layers.Dense(50, activation=tf.nn.relu),
    tf.keras.layers.Dropout(0.8),
    tf.keras.layers.Dense(1, activation=tf.nn.sigmoid)
])
model.compile(optimizer=moptimizer,loss='binary_crossentropy',metrics=['accuracy'])
history = model.fit(X_train, y_train, epochs=1)
    #load weights
model.load_weights(cp_filepath)


#retrain best with trva
trva_x=np.concatenate((X_train,X_valid),axis=0)
trva_y=np.concatenate((y_train,y_valid),axis=0)
tf.keras.backend.clear_session()
model = tf.keras.models.Sequential([
    tf.keras.layers.Dense(50, activation=tf.nn.relu),
    tf.keras.layers.Dropout(0.8),
    tf.keras.layers.Dense(50, activation=tf.nn.relu),
    tf.keras.layers.Dropout(0.8),
    tf.keras.layers.Dense(1, activation=tf.nn.sigmoid)
])
model.compile(optimizer=moptimizer,loss='binary_crossentropy',metrics=['accuracy'])
history = model.fit(trva_x, trva_y, epochs=1)
model.load_weights(cp_filepath)


history2 = model.fit(trva_x, trva_y, epochs=300)

plt.figure(figsize=(10, 6))
def plot_graphs(history, string):
  plt.plot(history.history[string])
  plt.xlabel("Epochs")
  plt.ylabel(string)
  plt.show()
plot_graphs(history2, "accuracy")
#plt.figure(figsize=(10, 6))#separate
plot_graphs(history2, "loss")
plt.legend(["accuracy", "loss"])

results = model.evaluate(x_test, y_test, batch_size=128)#0.74!