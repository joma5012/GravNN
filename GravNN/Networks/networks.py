import tensorflow as tf

def TraditionalNet(layers, activation, **kwargs):
    inputs = tf.keras.Input(shape=(layers[0],))
    x = inputs
    for i in range(1,len(layers)-1):
        x = tf.keras.layers.Dense(units=layers[i], 
                                    activation=activation, 
                                    kernel_initializer='glorot_normal',
                                    )(x)
                                    #dtype=kwargs['dtype'])(x)
        if 'dropout' in kwargs:
            if kwargs['dropout'] != 0.0:
                x = tf.keras.layers.Dropout(kwargs['dropout'])(x)
    outputs = tf.keras.layers.Dense(units=layers[-1], 
                                    activation='linear', 
                                    kernel_initializer='glorot_normal',
                                    dtype='float32'
                                    )(x)
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    return model

def ResNet(layers, activation, **kwargs):
    skip_offset = 3
    inputs = tf.keras.Input(shape=(layers[0],))
    x = inputs
    for i in range(1,len(layers)-1):
        x = tf.keras.layers.Dense(units=layers[i], 
                                    activation=activation, 
                                    kernel_initializer='glorot_normal')(x)
        if (i-1) % skip_offset == 0 and (i-1) == 0:
            skip = x 
        if (i-1) % skip_offset == 0 and (i-1) != 0:
            x = tf.keras.layers.Add()([x, skip])
            x = tf.keras.layers.Activation(activation)(x)
            skip = x 
    outputs = tf.keras.layers.Dense(units=layers[-1], 
                                    activation='linear', 
                                    kernel_initializer='glorot_normal')(x)
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    return model

def InceptionNet(layers, activation, **kwargs):
    inputs = tf.keras.Input(shape=(layers[0],))
    x = inputs
    for i in range(1,len(layers)-1):
        tensors = []
        for j in range(0, len(layers[i])):
            x_j = tf.keras.layers.Dense(units=layers[i][j], 
                                        activation=activation, 
                                        kernel_initializer='glorot_normal')(x)
            tensors.append(x_j)
        x = tf.keras.layers.Concatenate(axis=1)(tensors)
        x = tf.keras.layers.Activation(activation)(x)
        
    outputs = tf.keras.layers.Dense(units=layers[-1], 
                                    activation='linear', 
                                    kernel_initializer='glorot_normal')(x)
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    return model

def DenseNet(layers, activation, **kwargs):
    inputs = tf.keras.Input(shape=(layers[0],))
    x = inputs
    for i in range(1,len(layers)-1):
        tensors = []
        tensors.append(x)
        if len(layers[i]) > 1:
            for j in range(0, len(layers[i])):
                y = tf.keras.layers.Dense(units=layers[i][j], 
                                            activation=activation, 
                                            kernel_initializer='glorot_normal')(x)
                tensors.append(y)
                x = tf.keras.layers.Concatenate(axis=1)(tensors)
        else:
            x = tf.keras.layers.Dense(units=layers[i][0],
                                        activation=activation,
                                        kernel_initializer='glorot_normal')(x)

    outputs = tf.keras.layers.Dense(units=layers[-1], 
                                    activation='linear', 
                                    kernel_initializer='glorot_normal')(x)
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    return model