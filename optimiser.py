import numpy as np
import matplotlib.pyplot as plt

# Constants
g0 = 9.80665 # standard gravity (m/s^2)
sgp = 3.986*10**14 # standard gravitational parameter (m^3/s^2)
tstep = 0.1

# Functions
def calc_C3(position, velocity, sgp): #sgp = standard gravitational parameter, inputs/outputs in kgms
    return np.linalg.norm(velocity)**2 - 2*sgp/np.linalg.norm(position)

def vel_from_C3(C3, sgp, position):
    return np.sqrt(C3 + 2*sgp/np.linalg.norm(position))

def calc_grav(position, sgp): # returns acceleration due to gravity (m/s^2)
    return -sgp/np.linalg.norm(position)**2

# def calc_pay_frac(delta_v, isp, pmf):
#     return (np.exp(delta_v/(isp*g0))*(pmf - 1) + 1)/(np.exp(delta_v/(isp*g0))*pmf)

# Sim variables
position_0 = np.array([0 , 6771000]) # initial xy position (m)
velocity_0 = np.array([np.sqrt(sgp/np.linalg.norm(position_0)), 0]) # inital xy velocity (m/s)
traj_C3 = -1.67*10**6
traj_C3 = 8.19*10**6

class vehicle:

    def __init__(self, isp, eng_twr, pmf):
        self.isp = isp
        self.eng_twr = eng_twr
        self.pmf = pmf

    def __str__(self):
        return f"Isp(s): {self.isp} Engine TWR: {self.eng_twr} Propellant Mass Fraction: {self.pmf}"
    
    def sim(self, a0, pos, vel, target_C3):
        if a0/g0 > self.eng_twr:
            return 'ERROR: a0 higher than engine TWR'
        else:
            t = 0
            while calc_C3(pos, vel, sgp) < target_C3:
                a = a0/(1 - (a0/(self.isp*g0))*t) # instantaneous acceleration due to thrust (m/s^2)
                acceleration = a*vel/np.linalg.norm(vel) + calc_grav(pos, sgp)*pos/np.linalg.norm(pos)
                pos = pos + (vel + 0.5*acceleration)*tstep
                vel = vel + acceleration*tstep
                t += tstep
            return self.isp*g0*np.log(1/(1 - a0/(self.isp*g0)*t)) # returns expended Delta-V (m/s)
        
    def calc_pay_frac(self, delta_v, acceleration):
        payload = (np.exp(delta_v/(self.isp*g0))*(self.pmf - 1) + 1)/(np.exp(delta_v/(self.isp*g0))*self.pmf) - acceleration/(self.eng_twr*g0)
        dry = acceleration/(self.eng_twr*g0) + (1 - np.exp(-delta_v/(self.isp*g0)))*(1/self.pmf - 1)
        return (payload)
    
    def plot(self, pos, vel, target_C3):
        resolution = 50
        accelerations = np.linspace(0.5, 3, resolution)
        payload_fractions = np.zeros(resolution)
        delta_vs = np.zeros(resolution)
        drys = np.zeros(resolution)
        print("")
        print("running")
        print("")
        for i in range(resolution):
            delta_vs[i] = self.sim(accelerations[i], pos, vel, target_C3)
            payload_fractions[i] = self.calc_pay_frac(delta_vs[i]) - accelerations[i]/(self.eng_twr*g0) 

        plt.figure()
        plt.plot(accelerations, payload_fractions)
        plt.xlabel('Initial Acceleration (m/s^2)')
        plt.ylabel('Payload Fraction to C3 = ' + str(target_C3*10**-6) + 'km^2/s^2')
        plt.title('Isp = ' + str(self.isp) + 's, TWR = ' + str(self.eng_twr) + ', Pmf = ' + str(self.pmf))
        plt.grid()

        plt.figure()
        plt.plot(accelerations, delta_vs)
        plt.xlabel('Initial Acceleration (m/s^2)')
        plt.ylabel('Delta-V (m/s) to reach C3 = ' + str(target_C3*10**-6) + 'km^2/s^2')
        plt.title('Isp = ' + str(self.isp) + 's, TWR = ' + str(self.eng_twr) + ', Pmf = ' + str(self.pmf))
        plt.axhline(vel_from_C3(traj_C3, sgp, position_0) - np.linalg.norm(velocity_0), color = 'k', label = "Instantaneous Delta-V")
        plt.legend()
        plt.grid()

        index_max = payload_fractions.argmax()
        print(self)
        print(np.max(payload_fractions), accelerations[index_max], delta_vs[index_max], drys[index_max])
        print("")
        return (np.max(payload_fractions), accelerations[index_max], delta_vs[index_max], drys[index_max])

case_1 = vehicle(900, 2.14, 0.814)
case_2 = vehicle(460, 47.9, 0.897)
case_3 = vehicle(380, 250/1.5, 29/30)
case_1.plot(position_0, velocity_0, traj_C3)
#case_2.plot(position_0, velocity_0, traj_C3)
#case_3.plot(position_0, velocity_0, traj_C3)
plt.show()
