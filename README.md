# SDR_Matlab_OFDM_802.11n
[![npm](https://img.shields.io/npm/l/express.svg)](https://github.com/MeowLucian/SDR_Matlab_OFDM_802.11n/blob/master/LICENSE)

This simple OFDM demo is based on IEEE 802.11n OFDM format for testing SDR hardware.

Using Software Designed Radio (SDR) to transmit OFDM signals at 5 GHz.

Transmitter and Receiver hardware : Zedboard (Xilinx Zynq®-7000) + AD9361 (Analog Device-FMCOMMS3)

![Hardware](https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/Hardware.jpg)

# Software needed :

* Matlab R2016a
* support package for Xilinx Zynq-Based Radio (R2016a version)

Before executing this code, please follow this [tutorial](https://www.mathworks.com/help/supportpkg/xilinxzynqbasedradio/ug/install-support-package-for-xilinx-zynq-based-radio.html) to install support package for Xilinx Zynq-Based Radio in Matlab.

# Code Structure :
Please open multiple Matlab windows to run `Hardware_TX.m` and `Hardware_RX.m` respectively.

## Hardware_TX.m
> TX_signal.mat

> OFDM_TX.m
* data_Payload_1.mat
* data_Payload_2.mat
* oversamp.m
* setstate0_TX.m

## Hardware_RX.m
> OFDM_RX.m
* Long_preamble_slot_Frequency.mat
* HTL_k_slot_Frequency.mat
* setstate0_RX.m

## RX_test
* RX.mat
* RX2.mat

# GUI :
* GUI_TX

![Program GUI_TX](https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/GUI_TX.png)

* GUI_RX

![Program GUI_RX](https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/GUI_RX.png)

Video Demo : https://www.youtube.com/watch?v=CiHpd1Z8qK8
![Program GUI gif](https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/GUI_gif.gif)

# System Model :

## OFDM Block Diagram

<img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/OFDM_Block_Diagram.png" width="800">

## Code Function :

### Implemented

* Data Signal Mapping
* Packet Detection
* Coarse/Fine Frequency Offset Estimation & Compensation
* Channel Estimation & One-Tap Equalizer
* Data De-Mapping

### Not implemented yet

* AGC (Auto Gain Control)
* Fine Symbol Timing Estimation
* Convolutional Decoding
* De-Interleaving
* Frame check calculation (It assumes every packet is the same length)

## TX System Model
<img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/TX%20System%20Model.png" width="500">

> * Short Preamble
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/TX%20System%20Model_Short%20Preamble.png" width="800">

> * Long Preamble
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/TX%20System%20Model_Long%20Preamble.png" width="800">

> * HT-Long Preamble (Tone Interleaving)
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/TX%20System%20Model_HTLong%20Preamble.png" width="800">

> * Payload
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/TX%20System%20Model_Payload.png" width="600">

> * TX signal
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/TX%20System%20Model_Final.png" width="700">

## TX RX Hardware Parameters
| Center Frequency                 | 5 GHz                            |
|:--------------------------------:|:--------------------------------:|
| Baseband Sample Rate (Bandwidth) | 20 MHz                           |
| Ts (Sampling time)               | 50 ns                            |
| Samples Per Frame                | 3500                             |
| PC Host IP address               | 192.168.3.1                      |
| TX IP address                    | 192.168.3.2                      |
| RX IP address                    | 192.168.3.3                      |

## The way to change Hardware IP / Mac address

Edit `newip.sh` file in SD card

```
# Flush existing config
ip addr flush dev eth0
ip link set dev eth0 down
# Set up new config
ip addr add 192.168.3.3/24 dev eth0
ip link set eth0 address 00:0A:35:00:01:23
ip route add default via 192.168.3.1
ip link set dev eth0 up
```

Then, use router DHCP hand setting mode to distribute network configuration parameters :

![Router setting](https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/Router%20setting.PNG)

## RX System Model
<img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/RX%20System%20Model.png" width="300">

> * "Delay and Correlate" algorithm for Packet Detection
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/RX%20System%20Model_Delay%20and%20Correlate%20algorithm.png" width="500">

> * Packet Detection (normal case) , Threshold=0.75
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/RX%20System%20Model_Packet%20Detection.png" width="600">

> * Coarse CFO Estimation & Compensation
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/RX%20System%20Model_Coarse%20CFO%20Estimation.png" width="600">

> * Fine CFO Estimation & Compensation
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/RX%20System%20Model_Fine%20CFO%20Estimation.png" width="600">

> * MIMO Channel Estimation Diagram (Tone Interleaving)
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/RX%20System%20Model_Channel%20Estimation.png" width="600">

> * Channel Estimation & Equalizer
> <img src="https://raw.githubusercontent.com/MeowLucian/SDR_Matlab_OFDM_802.11n/master/Picture/RX%20System%20Model_Channel%20Estimation%20%26%20Equalizer.png" width="600">
