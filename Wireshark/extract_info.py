import numpy as np

def initialize(frame_type, frame_len, data_rate, duration, length, drate, frduration):
    length[frame_type][0] = frame_len
    length[frame_type][1] = frame_len

    drate[frame_type][0] = data_rate
    drate[frame_type][1] = data_rate

    frduration[frame_type][0] = duration
    frduration[frame_type][1] = duration
    return

def frame_type_resolved(frame_type):
    switcher = {
        0: "Association Request",
        1: "Association Response",
        2: "Reassociation Request",
        3: "Reassociation Response",
        4: "Probe Request",
        5: "Probe Response",
        6: "6",
        7: "7",
        8: "Beacon",
        9: "ATIM",
        10: "Disassociation",
        11: "Authentication",
        12: "Deauthentication",
        13: "Action Frames",
        14: "14",
        15: "15",
        16: "16",
        17: "17",
        18: "18",
        19: "19",
        20: "20",
        21: "21",
        22: "22",
        23: "23",
        24: "Block ACK Request",
        25: "Block ACK",
        26: "Power Save Poll",
        27: "Request to Send",
        28: "Clear to Send",
        29: "ACK",
        30: "CFP End",
        31: "CFP End ACK",
        32: "Data",
        33: "Data + CF ACK",
        34: "Data + CF Poll",
        35: "Data + CF ACK + CF Poll",
        36: "Null Data",
        37: "Null Data + CF ACK",
        38: "Null Data + CF Poll",
        39: "Null Data + CF ACK + CF Poll",
        40: "QoS Data",
        41: "QoS Data + CF ACK",
        42: "QoS Data + CF Poll",
        43: "QoS Data + CF ACK + CF Poll",
        44: "Null QoS Data",
        45: "45",
        46: "Null QoS Data + CF Poll",
        47: "Null QoS Data + CF ACK + CF Poll"
    }
    return switcher.get(frame_type, "nothing")

counterframes = [0] * 48
length = np.zeros((48, 4))
drate = np.zeros((48, 4))
frduration = np.zeros((48, 4))
fcsstatus = np.zeros((48, 2))
loop = 0
initial = [1] * 48
with open ('library.txt') as fo:
    for rec in fo:
        frame_type, frame_len, data_rate, duration, status = rec.split()

        if(loop == 1):
            frame_type = int(frame_type)
            frame_len = int(frame_len)
            data_rate = float(data_rate)
            duration = int(duration)
            status = int(status)
            
            counterframes[frame_type] += 1
            
            if(initial[frame_type] == 1):
                initialize(frame_type, frame_len, data_rate, duration, length, drate, frduration)
                initial[frame_type] = 0
            #column 0 is min, 1 is max, 2 is sum
            min = length[frame_type][0]
            max = length[frame_type][1]
            if(frame_len < min):
                length[frame_type][0] = frame_len
            if(frame_len > max):
                length[frame_type][1] = frame_len
            length[frame_type][2] = length[frame_type][2] + frame_len

            min = drate[frame_type][0]
            max = drate[frame_type][1]
            if(data_rate < min):
                drate[frame_type][0] = data_rate
            if(data_rate > max):
                drate[frame_type][1] = data_rate
            drate[frame_type][2] = drate[frame_type][2] + data_rate

            min = frduration[frame_type][0]
            max = frduration[frame_type][1]
            if(duration < min):
                frduration[frame_type][0] = duration
            if(duration > max):
                frduration[frame_type][1] = duration
            frduration[frame_type][2] = frduration[frame_type][2] + duration

            #column 0 is sum of good fcs, 1 is sum of bad fcs
            if(status == 1):
                fcsstatus[frame_type][0] += 1
            else:
                fcsstatus[frame_type][1] += 1
        loop = 1

    #column 3 is average
    for i in range(0, 48):
        if(counterframes[i] > 0):
            length[i][3] = length[i][2] / counterframes[i]
            drate[i][3] = drate[i][2] / counterframes[i]
            frduration[i][3] = frduration[i][2] / counterframes[i]

with open('library_nums.txt', 'w') as outfile:
    print('Frame_Type, Number, Min_Len, Max_Len, Av_Len, Min_drate, Max_drate, Av_drate, Min_dur, Max_dur, Av_dur, Good, Bad', file=outfile) 
    for i in range(0, 48):
        if(counterframes[i] > 0):
            print(i, '\t', counterframes[i], '\t', int(length[i][0]), '\t', int(length[i][1]), '\t', format(length[i][3], '.2f'), '\t', \
                  format(drate[i][0], '.2f'), '\t', format(drate[i][1], '.2f'), '\t', format(drate[i][3], '.2f'), '\t', \
                  int(frduration[i][0]), '\t', int(frduration[i][1]), '\t', format(frduration[i][3], '.2f'), '\t', int(fcsstatus[i][0]), '\t', int(fcsstatus[i][1]), file=outfile)

with open('library_names.txt', 'w') as outfile:
    print('Frame_Type, Number, Min_Len, Max_Len, Av_Len, Min_drate, Max_drate, Av_drate, Min_dur, Max_dur, Av_dur, Good, Bad', file=outfile) 
    for i in range(0, 48):
        if(counterframes[i] > 0):
            fr_type_res = frame_type_resolved(i)
            print(fr_type_res, '\t', counterframes[i], '\t', int(length[i][0]), '\t', int(length[i][1]), '\t', format(length[i][3], '.2f'), '\t', \
                  format(drate[i][0], '.2f'), '\t', format(drate[i][1], '.2f'), '\t', format(drate[i][3], '.2f'), '\t', \
                  int(frduration[i][0]), '\t', int(frduration[i][1]), '\t', format(frduration[i][3], '.2f'), '\t', int(fcsstatus[i][0]), '\t', int(fcsstatus[i][1]), file=outfile)
            
