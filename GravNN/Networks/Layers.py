from asyncio import constants
import tensorflow as tf
import numpy as np


class PreprocessingLayer(tf.keras.layers.Layer):
    def __init__(self, min, scale, dtype):
        super(PreprocessingLayer, self).__init__(dtype=dtype)
        self.min = tf.constant(min, dtype=dtype)
        self.scale = tf.constant(scale, dtype=dtype)

    def call(self, inputs):
        normalized_inputs = inputs * self.scale + self.min
        return normalized_inputs

    def get_config(self):
        config = super().get_config().copy()
        config.update(
            {
                "min": self.min,
                "scale": self.scale,
            }
        )
        return config

class PostprocessingLayer(tf.keras.layers.Layer):
    def __init__(self, min, scale, dtype):
        super(PostprocessingLayer, self).__init__(dtype=dtype)
        self.min = tf.constant(min, dtype=dtype)
        self.scale = tf.constant(scale, dtype=dtype)

    def call(self, normalized_inputs):
        inputs = (normalized_inputs - self.min)/self.scale
        return inputs

    def get_config(self):
        config = super().get_config().copy()
        config.update(
            {
                "min": self.min,
                "scale": self.scale,
            }
        )
        return config

class Cart2PinesSphLayer(tf.keras.layers.Layer):
    def __init__(self, dtype):
        """Successor to the Cart2SphLayer. The layer takes a
        cartesian input and transforms it into a non-singular spherical
        representation (see Pines formulation). This bypasses the singularity introduced at the pole
        when taking a derivative of the potential.

        https://ntrs.nasa.gov/api/citations/19760011100/downloads/19760011100.pdf
        defines of alpha, beta, and gamma (i.e. three angle non-singular system)
        """
        super(Cart2PinesSphLayer, self).__init__(dtype=dtype)

    def call(self, inputs):
        inputs_transpose = tf.transpose(inputs)
        X = inputs_transpose[0]
        Y = inputs_transpose[1]
        Z = inputs_transpose[2]

        XX = tf.square(X)
        YY = tf.square(Y)
        ZZ = tf.square(Z)
        r = tf.sqrt(XX + YY + ZZ)  # r

        s = X / r  # sin(beta)
        t = Y / r  # sin(gamma)
        u = Z / r  # sin(alpha)

        spheres = tf.stack([r, s, t, u], axis=1)
        return spheres

    def get_config(self):
        config = super().get_config().copy()
        return config

class PinesSph2NetLayer(tf.keras.layers.Layer):
    def __init__(self, dtype, scalers, mins, ref_radius):
        """Layer used to normalize the r value of the Cart2PinesSphLayer
        between the values of [-1,1]

        Args:
                scalers (np.array): values used to scale each input quantity. These are computed before initialization
                mins (np.array): values used to bias each input before scaling. These are computed before initialization
        """
        super(PinesSph2NetLayer, self).__init__(dtype=dtype)
        self.scalers = tf.constant(scalers, dtype=dtype).numpy()
        self.mins = tf.constant(mins, dtype=dtype).numpy()
        self.ref_radius = tf.constant(ref_radius, dtype=dtype).numpy()

    def call(self, inputs):
        inputs_transpose = tf.transpose(inputs)
        r = inputs_transpose[0]
        beta = inputs_transpose[1]
        gamma = inputs_transpose[2]
        alpha = inputs_transpose[3]

        # Normalize r -> [-1,1]
        r = tf.add(tf.multiply(r, self.scalers[0]), self.mins[0])
        spheres = tf.stack([r, beta, gamma, alpha], axis=1)
        return spheres

    def get_config(self):
        config = super().get_config().copy()

        config.update(
            {
                "scalers": self.scalers,
                "mins": self.mins,
                "ref_radius":self.ref_radius
            }
        )
        config = super().get_config().copy()
        return config

class PinesSph2NetRefLayer(tf.keras.layers.Layer):
    def __init__(self, dtype, scalers, mins, ref_radius):
        """Layer used to normalize the r value of the Cart2PinesSphLayer
        between the values of [-1,1]

        Args:
                scalers (np.array): values used to scale each input quantity. These are computed before initialization
                mins (np.array): values used to bias each input before scaling. These are computed before initialization
        """
        super(PinesSph2NetRefLayer, self).__init__(dtype=dtype)
        self.ref_radius = ref_radius

    def call(self, inputs):
        inputs_transpose = tf.transpose(inputs)
        r = inputs_transpose[0]
        beta = inputs_transpose[1]
        gamma = inputs_transpose[2]
        alpha = inputs_transpose[3]

        # Normalize r -> [0,inf]
        r_inv_prime = tf.divide(self.ref_radius,r)
        spheres = tf.stack([r_inv_prime, beta, gamma, alpha], axis=1)
        return spheres

    def get_config(self):
        config = super().get_config().copy()

        config.update(
            {
                "ref_radius":self.ref_radius
            }
        )
        config = super().get_config().copy()
        return config

class PinesSph2NetLayer_v2(tf.keras.layers.Layer):
    def __init__(self, dtype, ref_radius_min=None, ref_radius_max=None):
        super(PinesSph2NetLayer_v2, self).__init__(dtype=dtype)
        self.ref_radius_min = tf.cond(ref_radius_min != None, 
                        lambda : tf.constant(ref_radius_min, dtype=dtype), 
                        lambda : tf.constant(1E-3, dtype=dtype)).numpy() 
                        # Don't set min to zero otherwise nans
        self.ref_radius_max = tf.cond(ref_radius_max != None, 
                        lambda : tf.constant(ref_radius_max, dtype=dtype), 
                        lambda : tf.constant(1.0, dtype=dtype)).numpy()
        
        # bounds of the feature -- shouldn't necessarily be a large difference.
        # think about how much the output should change with respect to the inputs
        self.s_min = tf.constant(1.0, dtype=dtype).numpy()
        self.s_max = tf.constant(1.0 + (self.ref_radius_max - self.ref_radius_min)/self.ref_radius_min, dtype=dtype).numpy()           

    def call(self, inputs):
        inputs_transpose = tf.transpose(inputs)
        r = inputs_transpose[0]
        s = inputs_transpose[1]
        t = inputs_transpose[2]
        u = inputs_transpose[3]

        scale = tf.divide(self.s_max - self.s_min, self.ref_radius_max - self.ref_radius_min)
        min_arg = self.s_min - tf.multiply(self.ref_radius_min, scale)
        r_star = tf.multiply(r, scale) + min_arg

        one = tf.constant(1.0, dtype=r.dtype)
        r_inv = tf.divide(one, r_star)
        # r_inv = r_star
        spheres = tf.stack([r_inv, s, t, u], axis=1)
        return spheres

    def get_config(self):
        config = super().get_config().copy()
        config.update(
            {
                "dtype" : self.dtype,
                "s_min" : self.s_min,
                "s_max" : self.s_max,
                "ref_radius_min" : self.ref_radius_min,
                "ref_radius_max" : self.ref_radius_max
            }
        )
        return config


class AugmentedPotentialLayer(tf.keras.layers.Layer):
    def __init__(self, dtype, mu, r_max):
        super(AugmentedPotentialLayer, self).__init__(dtype=dtype)
        self.mu = tf.constant(mu, dtype=dtype).numpy()
        self.r_max = tf.constant(r_max, dtype=dtype).numpy()

    def call(self, u_nn, inputs):
        one = tf.constant(1.0, dtype=u_nn.dtype)
        half = tf.constant(0.5, dtype=u_nn.dtype)
        k = tf.constant(10000.0, dtype=u_nn.dtype)
        r = tf.linalg.norm(inputs, axis=1, keepdims=True)
        dr = tf.subtract(r,self.r_max)
        h = half+half*tf.tanh(k*dr)
        u_pm = tf.negative(tf.divide(self.mu,r))
        u_model = (one - h)*(u_nn + u_pm) + h*u_pm 
        return u_model

    def get_config(self):
        config = super().get_config().copy()
        config.update(
            {
                "mu" : self.mu,
                "r_max" : self.r_max,
            }
        )
        return config


class PointMassLayer(tf.keras.layers.Layer):
    def __init__(self, dtype, mu, r_max):
        super(PointMassLayer, self).__init__(dtype=dtype)
        self.mu = tf.constant(mu, dtype=dtype).numpy()

    def call(self, u_nn, inputs):
        r = tf.linalg.norm(inputs, axis=1, keepdims=True)
        u_pm = tf.negative(tf.divide(self.mu,r))
        return u_pm

    def get_config(self):
        config = super().get_config().copy()
        config.update(
            {
                "mu" : self.mu,
            }
        )
        return config


class BlendPotentialLayer(tf.keras.layers.Layer):
    def __init__(self, dtype, mu, r_max):
        super(BlendPotentialLayer, self).__init__(dtype=dtype)
        self.mu = tf.constant(mu, dtype=dtype).numpy()
        self.r_max = tf.constant(r_max, dtype=dtype).numpy()

    def build(self, input_shapes):
        self.radius = self.add_weight("radius",
                            shape=[1],
                            trainable=True, 
                            initializer =tf.keras.initializers.Constant(value=self.r_max),
                            )
        self.k = self.add_weight("k",
                            shape=[1],
                            trainable=True, 
                            initializer =tf.keras.initializers.Constant(value=10),
                            )
        super(BlendPotentialLayer, self).build(input_shapes)


    def call(self, u_nn, u_analytic, inputs):
        one = tf.constant(1.0, dtype=u_nn.dtype)
        half = tf.constant(0.5, dtype=u_nn.dtype)
        r = tf.linalg.norm(inputs, axis=1, keepdims=True)
        dr = tf.subtract(r,self.radius)
        h = half+half*tf.tanh(self.k*dr)
        u_model = (one - h)*(u_nn + u_analytic) + h*u_analytic 
        return u_model

    def get_config(self):
        config = super().get_config().copy()
        config.update(
            {
                "mu" : self.mu,
                "r_max" : self.r_max,
            }
        )
        return config




class PinesAlgorithmLayer(tf.keras.layers.Layer):
    def __init__(self, dtype, mu, a, cBar, sBar):
        super(PinesAlgorithmLayer, self).__init__(dtype=dtype)
        self.mu = tf.constant(mu, dtype=dtype).numpy()
        self.a = tf.constant(a, dtype=dtype).numpy()
        self.cBar = tf.constant(cBar, dtype=dtype).numpy()
        self.sBar = tf.constant(sBar, dtype=dtype).numpy()
        self.N = tf.constant(len(cBar)-3, dtype=tf.int32).numpy()
        self.n1, self.n2 = self.compute_normalization_constants(self.N)
        # a = self.compute_aBar(tf.constant(10,dtype=dtype))
        # rE, iM = self.compute_rE_iM(tf.constant(10,dtype=dtype), tf.constant(10,dtype=dtype))

    def getK(self, x):
        return tf.constant(1.0, dtype=self.dtype) if (x == 0) else tf.constant(2.0, dtype=self.dtype)

    def compute_normalization_constants(self, N):
        n1 = tf.zeros((N + 2, N + 2), dtype=self.dtype) 
        n2 = tf.zeros((N + 2, N + 2), dtype=self.dtype) 

        for l_idx in range(0, N + 2):
            for m in range(0, l_idx + 1):
                if l_idx >= m + 2:
                    l = tf.constant([l_idx], self.dtype)
                    n1_lm = tf.sqrt(
                        ((2.0 * l + 1.0) * (2.0 * l - 1.0)) / ((l - m) * (l + m))
                    )
                    n2_lm = tf.sqrt(
                        ((l + m - 1.0) * (2.0 * l + 1.0) * (l - m - 1.0))
                        / ((l + m) * (l - m) * (2.0 * l - 3.0))
                    )
                    n1 = tf.tensor_scatter_nd_update(n1, [[l_idx,m]], n1_lm, name='n1_update')
                    n2 = tf.tensor_scatter_nd_update(n2, [[l_idx,m]], n2_lm, name='n2_update')

        return n1.numpy(), n2.numpy()

    def compute_rE_iM(self, s, t):
        rE = tf.scatter_nd(tf.constant([[0]]), tf.constant([1.0], dtype=self.dtype), shape=tf.constant([self.N+2]), name='rE')
        iM = tf.scatter_nd(tf.constant([[0]]), tf.constant([0.0], dtype=self.dtype), shape=tf.constant([self.N+2]), name='iM')
        
        for i in range(1, self.N+2):
            rE_m1 = tf.gather(rE, [i-1])
            iM_m1 = tf.gather(iM, [i-1])
            rE = tf.tensor_scatter_nd_update(rE, [[i]], s * rE_m1 - t * iM_m1, name='rE_update') # introduces error
            iM = tf.tensor_scatter_nd_update(iM, [[i]], s * iM_m1 + t * rE_m1, name='iM_update')

        return rE, iM

    def compute_aBar(self, u):
        N = self.N
        aBar = tf.scatter_nd([[0,0]], tf.constant([1.0], dtype=u.dtype), shape=((N+2, N+2)), name='aBar')

        for l in tf.range(1, N + 2):
            a_lm1_lm1 = tf.gather_nd(aBar, [[l-1,l-1]], name='aBar_gather')
            l_float = tf.cast(l, dtype=self.dtype)
            a_l_l = tf.sqrt((2.0 * l_float + 1.0) * self.getK(l) / (2.0 * l_float * self.getK(l - 1))) * a_lm1_lm1
            a_l_lm1 = tf.sqrt(2.0 * l_float * self.getK(l - 1) / self.getK(l)) * a_l_l * u
            aBar = tf.tensor_scatter_nd_update(aBar, [[l,l]], a_l_l, name='aBar_update_1')
            aBar = tf.tensor_scatter_nd_update(aBar, [[l,l-1]], a_l_lm1, name="aBar_update_2")

        for m in range(0, N + 2):
            for l in range(m + 2, N + 2):
                a_lm1_m = tf.gather_nd(aBar, [[l-1, m]], name='a_lm1_m')
                a_lm2_m = tf.gather_nd(aBar, [[l-2, m]], name='a_lm2_m')
                n1_lm = tf.gather_nd(self.n1, [[l,m]], name='n1_lm')
                n2_lm = tf.gather_nd(self.n2, [[l,m]], name='n2_lm')
                a_lm = u * n1_lm * a_lm1_m - n2_lm * a_lm2_m
                aBar = tf.tensor_scatter_nd_update(aBar, [[l,m]], a_lm, name='aBar_final')

        return aBar 

    def compute_rhol(self, a, r):
        rho = a / r
        rhol = tf.zeros((self.N+1), dtype=self.dtype)
        rhol = tf.scatter_nd([[0]], [self.mu/r], shape=((self.N+1,)), name='rho')
        rhol = tf.tensor_scatter_nd_update(rhol, [[0]], [self.mu/r], name='rho_update_0') # good
        rhol = tf.tensor_scatter_nd_update(rhol, [[1]], [self.mu/r * rho], name='rho_update_1') # good
        for l in range(1, self.N):
            rho_i = tf.gather(rhol,[l], name='rho_gather')
            rhol = tf.tensor_scatter_nd_update(rhol, [[l + 1]], rho * rho_i, name='rho_update_2') # introduce error
        return rhol

    def compute_potential(self, r, s, t, u, a):
        rE, iM = self.compute_rE_iM(s,t)
        rhol = self.compute_rhol(a, r)
        aBar = self.compute_aBar(u)

        potential = 0.0
        for l in range(1, self.N + 1):
            for m in range(0, l + 1):
                potential += rhol[l] * aBar[l,m] * \
                    (self.cBar[l,m] * rE[m] + self.sBar[l,m] * iM[m])

        potential += self.mu / r
        neg_potential = tf.negative(potential)
        return neg_potential


    def call(self, inputs):
        inputs_transpose = tf.transpose(inputs)
        r = inputs_transpose[0]
        s = inputs_transpose[1]
        t = inputs_transpose[2]
        u = inputs_transpose[3]
        a = tf.ones_like(r)*self.a

        # potential = self.compute_potential(r, s, t, u, a)

        potential = tf.map_fn(
            lambda x: self.compute_potential(
                x[0], x[1], x[2], x[3], x[4]), 
            elems=(r, s, t, u, a),
            fn_output_signature=(r.dtype),
            # parallel_iterations=10
            )
        u = tf.reshape(potential, (-1,1))
        return u


    def get_config(self):
        config = super().get_config().copy()
        config.update(
            {
                "mu" : self.mu,
                "a" : self.a,
                "n1" : self.n1,
                "n2" :  self.n2,
                "cBar" : self.cBar,
                "sBar" : self.sBar,
                "N" : self.N,
                "dtype" : self.dtype,
            }
        )
        return config