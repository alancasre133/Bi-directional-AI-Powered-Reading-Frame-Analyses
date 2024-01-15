import pandas as pd
import numpy as np
import tensorflow as tf
from matplotlib import pyplot as plt
class Model:
    def __init__(self) -> None:
        pass
    def build_model(self,my_learning_rate):
        # Most simple tf.keras models are sequential. 
        # A sequential model contains one or more layers.
        model = tf.keras.models.Sequential()

        # Describe the topography of the model.
        # The topography of a simple linear regression model
        # is a single node in a single layer. 
        model.add(tf.keras.layers.Dense(units=1, 
                                        input_shape=(1,)))

        # Compile the model topography into code that 
        # TensorFlow can efficiently execute. Configure 
        # training to minimize the model's mean squared error. 
        model.compile(optimizer=tf.keras.optimizers.experimental.RMSprop(learning_rate=my_learning_rate),
                        loss="mean_squared_error",
                        metrics=[tf.keras.metrics.RootMeanSquaredError()])

        return model           


    def train_model(self,model, feature, label, epochs, batch_size):
        # Feed the feature values and the label values to the 
        # model. The model will train for the specified number 
        # of epochs, gradually learning how the feature values
        # relate to the label values. 
        history = model.fit(x=feature,
                            y=label,
                            batch_size=batch_size,
                            epochs=epochs)

        # Gather the trained model's weight and bias.
        trained_weight = model.get_weights()[0]
        trained_bias = model.get_weights()[1]

        # The list of epochs is stored separately from the 
        # rest of history.
        epochs = history.epoch
        
        # Gather the history (a snapshot) of each epoch.
        hist = pd.DataFrame(history.history)

        # Specifically gather the model's root mean 
        # squared error at each epoch. 
        rmse = hist["root_mean_squared_error"]

        return trained_weight, trained_bias, epochs, rmse
    def plot_the_model(self,trained_weight, trained_bias, feature, label):
        # Label the axes.
        plt.xlabel("feature")
        plt.ylabel("label")

        # Plot the feature values vs. label values.
        plt.scatter(feature, label)

        # Create a red line representing the model. The red line starts
        # at coordinates (x0, y0) and ends at coordinates (x1, y1).
        x0 = 0
        y0 = trained_bias
        x1 = feature[-1]
        y1 = trained_bias + (trained_weight * x1)
        plt.plot([x0, x1], [y0, y1], c='r')

        # Render the scatter plot and the red line.
        plt.show()

    def plot_the_loss_curve(self,epochs, rmse):
        plt.figure()
        plt.xlabel("Epoch")
        plt.ylabel("Root Mean Squared Error")

        plt.plot(epochs, rmse, label="Loss")
        plt.legend()
        plt.ylim([rmse.min()*0.97, rmse.max()])
        plt.show()
class bitsequence:
        def __init__(self) ->None:
            self.length=0
            self.encode=0
        def _encode(self,nucl,i) -> int:
            base = 0
            if nucl == 'C' or nucl=='c':
                base = 10
            elif nucl == 'G' or nucl=='g':
                base = 20
            elif nucl == 'T' or nucl=='t':
                base = 30
            return 1 << (base + i)
        def fullencode(self,nucl) -> int:
            self.encode=0
            self.length=len(nucl)
            for i in range(len(nucl)):
                self.encode |=self._encode(nucl[i],i)
            return self.encode
        def decode(self,i):
            res = 'A'
            base = i // 10
            if base == 1:
                res = 'C'
            elif base == 2:
                res = 'G'
            elif base == 3:
                res = 'T'
            
            return res
        def fulldecode(self):
            decoded=['']*self.length
            count=0
            print(self.encode)
            for i in range(self.encode*400):
                if (1 << i) & self.encode:
                    decoded[i%self.length] = self.decode(i)
                       
                    count += 1
                if count == self.length:
                    break
            
            return ''.join(decoded)
from collections import defaultdict

code=["atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
      "atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
      "atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
      "atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
      "atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
      "atgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtga",
      "atgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtga",
      "atgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtga",
      "atgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtga",
      "atgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtga",
      "atgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtga",
      "atgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtga"
      ]
label=["atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
       "atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
       "atgacgtatccaaggaggcgttaccggagaagaagacaccgcccccgcagccatcttggccagatcctccgccgccgcccctggctcgtccacccccgccaccgttaccgctggagaaggaaaaatggcatcttcaacacccgcctctcccgcaccttcggatatactatcaagcgaaccacagtcaaaacgccctcctgggcggtggacatgatgagattcaatattaatgactttcttcccccaggagggggctcaaacccccgctctgtgccctttgaatactacagaataagaaaggttaaggttgaattctggccctgctccccgatcacccagggtgacaggggagtgggctccagtgctgttattctagatgataactttgtaacaaaggccacagccctcacctatgacccctatgtaaactactcctcccgccataccataacccagcccttctcctaccactcccgctactttacccccaaacctgtcctagattccactattgattacttccaaccaaacaacaaaagaaatcagctgtggctgagactacaaactactggaaatgtagaccacgtaggcctcggcactgcgttcgaaaacagtatatacgaccaggaatacaatatccgtgtaaccatgtatgtacaattcagagaatttaatcttaaagaccccccacttaacccctaa",
       "atggtaaccatcccaccacttgtttctaggtggtttccagtatgtggtttccgggtctgcaaaattagcagcccatttgcttttaccacacccaggtggccccacaatgacgtgtacattagtcttccaatcacgcttctgcattttcccgctcactttcaaaagttcagccagcccgcggaaatttctgacaaacgttacagggtgctgctctgcaacggtcaccagactcccgctctccaacaaggtactcacagcagtagacaggtcactccgttgtccctgagatctaggagctccacactccatcagtaa",
       "atgcccagcaagaagaatggaagaagcggaccccaaccccataaaaggtgggtgttcacactgaataatccttccgaagacgagcgcaagaaaatacgggatcttccaatatccctatttgattattttattgttggcgaggagggtaatgaggaaggacgaacacctcacctccaggggttcgctaattttgtgaagaagcagacttttaataaagtgaagtggtatttgggtgcccgctgccacatcgagaaagcgaaaggaacagatcagcagaataaagaatactgcagtaaagaaggcaacttactgatggagtgtggagctcctagatctcagggacaacggagtgacctgtctactgctgtgagtaccttgttggagagcgggagtctggtgaccgttgcagagcagcaccctgtaacgtttgtcagaaatttccgcgggctggctgaacttttgaaagtgagcgggaaaatgcagaagcgtgattggaagactaatgtacacgtcattgtggggccacctgggtgtggtaaaagcaaatgggctgctaattttgcagacccggaaaccacatactggaaaccacctagaaacaagtggtgggatggttaccatggtgaagaagtggttgttattgatgacttttatggctggctgccctgggatgatctactgagactgtgtgatcgatatccattgactgtagagactaaaggtggaactgtaccttttttggcccgcagtattctgattaccagcaatcagaccccgttggaatggtactcctcaactgctgtcccagctgtagaagctctttatcggaggattacttccttggtattttggaagaatgctacagaacaatccacggaggaagggggccagttcgtcaccctttcccccccatgccctgaatttccatatgaaataaattactga",
       "atgatatccatcccaccacttatttctactaggcttccagtaggtgtcgctaggctcagcaaaattacgggcccactggctcttcccacaaccgggcgggcccactatgacgtgtacagctgtcttccaatcacgctgctgcatcttcccgctcactttcaaaagttcagccagcccgcggaaatttctcacatacgttacagggaactgctcggctacagtcaccaaagaccccgtctccaaaagggtactcacagcagtagacaagtcgctgcgcttcccctggttccgcggagctccacactcgataagtatgtggccttctttactgcagtattctttattctgctggtcggttcctttcgctttctcgatgtggcagcgggcaccaaaataccacttcaccttgttaaaagtctgcttcttagcaaaattcgcaaacccctggaggtgaggagttctaccctcttccaaaccttcctctccgcaaacaaaataatcaaaaagggagattggaagctcccgtattttgtttttctcctcctcggaaggattattaagggtgaacacccacctcttatggggttgcgggccgcttttcttgcttggcattttcactga",
       "atgacgtggccaaggaggcgttaccgcagaagaaggacccgcccccgcagccatcttggaaacatcctccggagaagaccatatttggcacaccccgccttcagaaaccgttacagatggcgccgaaagacgggtatcttcaattcccgcctttctagagaatttgtactcaccataaaaggaggatactcgcagccatcttggaatgttaaccacctcaaattcaacatcggccagttcctccccccctcaggcggcaccaacccactacccctacctttccaatactaccgtattagaaaggctaaatatgaattttaccccagagaccccatcacctctaatcaaagaggtgttgggtccactgttgttatcttggatgccaactttgtaaccccctccaccaacttggcctatgacccctatattaactactcctcccgccacaccataaggcagccctttacctaccactccaggtacttcacccccaaacctgagctggaccaaacaattgattggttccacccaaataataaaagaaaccagctgtggctccatttaaatacccacaccaatgtcgagcacacaggcctcggctatgcgctccaaaatgcagccacagcccaaaattatgtggtaaggctgactatttatgtacaattcagagaatttatcctaaaagaccctctaaataaataa",
       "atgacaataaaagtccttatcttcaggacactcgtagcaccacaaaaactcagttatccagtccccgtcctgcggcatcaaaacacgaggccgcttcatcatccactgcccggtatatccactcacacggtgcggagcatccataatgggataccactttttctccctacagacctccgtggatccgggttccgaggtgccgggtaa",
       "atgcgagggcgtttacctgtgcccgcacccgaagcgcagcgggagcgcgcgcgaggggacacggcttgtcgccaccggaggggtcagatttatatttattttcacttagagaacggacttgtaacgaatccaaacttctttggtgccgtagaagtctgtcattccagttttttccgggacataaatgctccaaagcagtgctccccattgaacggtggggtcatatgtgttgagccatggggtgggtctggagaaaaagaagaggctttgtcctgggtgagcgctggtagttcccgccagaagtggtttgggggtgaagtaacggctgtgtttttttttagaagtcataactttacgagtggaactttccgcataagggtcgtcttggagccaagtgtttgtggtccaggcgccgtctag",
       "atgtttttattggtcctccaggatgcgggaaaacgcgggaagcttgtgcggatgcggctgcgcgggaattgcagttgtatttcaagccacgggggccttggtgggatggttataatggggagggtgctgttattttggatgatttttatgggtgggttccatttgatgaattgctga",
       "atgttgtaccggaggagtggtattcatcggagaatattcgtggaaagttggaggccttgtttaggaggttcactaaggttgtttgttggggggaggggggggtaa",
       "atgtccaccgcccaggagggcgttttgactgtggttcgcttgatagtatatccgaaggtgcgggagaggcgggtgttgaagatgccatttttccttctccagcggtaa"
       ]
index=0
feature=[]
labels=[]
f=bitsequence()
"""
for i in range(len(label)):
    feature.append(f.fullencode(code[i]))
    labels.append(f.fullencode(label[i]))
    index+=1
data=pd.DataFrame(np.array([feature,labels]))
data.corr()
print(feature)
print(labels)
"""
print(f"sequence: atgccaatttga encode: {float(f.fullencode('atgccaatttga'))}")

model=Model()
learning_rate=0.01
epochs=10
my_batch_size=4

my_model = model.build_model(learning_rate)
trained_weight, trained_bias, epochs, rmse = model.train_model(my_model, feature, 
                                                         labels, epochs,
                                                         my_batch_size)

model.plot_the_loss_curve(epochs, rmse)