from __future__ import absolute_import

import numpy as np
import struct as st
from datetime import datetime, timedelta
import pandas
import matplotlib.pyplot as plt

__all__ = ["TMFrame", "TMStream", "parse_tm_file"]

header_byte_length = 16
sas_header_byte_length = header_byte_length + 2  # HEROES header + the SAS sync byte
packet_byte_length = 110
data_byte_length = packet_byte_length - packet_byte_length

header_unpack_string = 'HBBHHII'
sas_header_unpack_string = header_unpack_string + 'H'
data_unpack_string = 'IBHHH'
#data_unpack_string = 'I'

sas1_sync_word = '0xeb90'
sas2_sync_word = '0xf626'

class TMFrame:
    def __init__(self, data):
        self.data = data
        self.time = None
        self.valid = False
        if self.is_valid():
            header = st.unpack(sas_header_unpack_string, data[0:sas_header_byte_length])
            parsed_data = st.unpack(data_unpack_string, data[sas_header_byte_length:sas_header_byte_length+12])
            self.sequence_number = parsed_data[0]
            self.time = datetime.fromtimestamp(header[6]) + timedelta(seconds=header[5]/1e9)
            self.image_max = data[95]
            self.valid = True
            self.housekeeping = [parsed_data[3], parsed_data[4]]
            self.ctl = st.unpack('ff', data[96:104])
    def is_valid(self):
        """Test if the data represents a valid packet"""
        if len(self.data) == packet_byte_length:
            header = st.unpack(sas_header_unpack_string, self.data[0:sas_header_byte_length]) 
            con2 = hex(header[0]) == '0xc39a'
            con3 = hex(header[1]) == '0x70'
            con4 = hex(header[2]) == '0x30'
            con5 = (hex(header[7]) == sas1_sync_word) or (hex(header[7]) == sas2_sync_word)
            return con2 and con3 and con4 and con5
        else:
            return False
        
class TMStream:
    def __init__(self, frames):
        self.frames = frames
        self.times = [frame.time for frame in frames]
    def frames_valid(self):
        if isinstance(self.frame[0], TMFrame([0])):
            return True
    def image_max(self):
        image_max = [frame.image_max for frame in self.frames]
        lc = pandas.Series(image_max, self.times)
        return lc
    def temperatures(self):
        result = []
        for i in np.arange(7):
            temp_sensor = [frame.housekeeping[0] for frame in self.frames if ((frame.sequence_number % 8) == i) ]
            temp_sensor_times = [frame.time for frame in self.frames if ((frame.sequence_number % 8) == i)]
            result.append(pandas.Series(temp_sensor, temp_sensor_times))
        return result

def parse_sastm_file(filename):
    data = np.fromfile(file, dtype=np.uint8)
    num_bytes_in_file = len(data)
    count_packets = 0
    frames = []
    index = 0
    prevFrame = TMFrame([0])    # dummy frame
    while(index < num_bytes_in_file - packet_byte_length-1):
            header = st.unpack(sas_header_unpack_string, data[index:index+sas_header_byte_length])
            if ((hex(header[0]) == '0xc39a') and (hex(header[1]) == '0x70') and (hex(header[2]) == '0x30') and (hex(header[7]) == sas1_sync_word)):
                currentFrame = TMFrame(data[index:index+packet_byte_length])
                if currentFrame.valid:
                    count_packets = count_packets + 1
                    index = index + packet_byte_length
                    if currentFrame.time != prevFrame.time:
                        frames.append(currentFrame)
                    prevFrame = currentFrame
            else:
                index = index + 1            
    return frames