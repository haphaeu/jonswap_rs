#
# Plot output from JONSWAP
#
# First run jonswap and output to out.txt
#
import pandas as pd
import matplotlib.pyplot as plt

with open('out.txt') as f:
    lines = f.readlines()
 
# Wave parameters
hs = lines[0].split(':')[-1].strip()
tp = lines[1].split(':')[-1].strip()
wy = lines[2].split(':')[-1].strip()
    
# Plot spectra
start = 4 + lines.index('Spectrum\n')
end = -1 + lines.index('Time History\n')
headers = lines[start - 2].split()
data = [list(map(float, l.split())) for l in lines[start:end]]
df = pd.DataFrame(data, columns=headers)
df.set_index('T', inplace=True)
df.loc[:, 'PM':'JS'].plot()
plt.xlabel('Period [s]')
plt.ylabel('Energy [m²/(rad/s)]')
plt.title(f'Wave Spectrum - Hs {hs}, Tp {tp}, gamma {wy}')
plt.legend()
plt.grid()
plt.savefig('spectra.png')

# Plot timetrace

headers = lines[end + 3].split()
data = [list(map(float, l.split())) for l in lines[end+4:]]
df = pd.DataFrame(data, columns=headers)
df.set_index('Time', inplace=True)
df.plot()
plt.xlabel('Time [s]')
plt.ylabel('Elevation [m)]')
plt.title(f'Wave Realisation - Hs {hs}, Tp {tp}, gamma {wy}')
plt.legend()
plt.grid()
plt.savefig('timetrace.png')
