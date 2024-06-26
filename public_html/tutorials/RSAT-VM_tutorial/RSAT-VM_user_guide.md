---
title: "RSAT Virtual machine (VM) User Guide"
author: "Jacques van Helden"
output:
  html_document:
    highlight: tango
    theme: united
    toc: yes
    toc_depth: 3
  pdf_document:
    highlight: zenburn
    toc: yes
    toc_depth: 3
css: ../../course.css
---

------------------------------------

## Introduction

This documents contains instructions to download, install, configure and use a stand-alone instance of the Regulatory Sequence Analysis Tools (<span class="concept">RSAT</span>) in the form of a Virtual Machine (<span class="concept">VM</span>) running under [VirtualBox](https://www.virtualbox.org/).

------------------------------------

## Requirements

The RSAT Virtual Machine (<span class="concept">RSAT-VM</span>) requires

1. The [VirtualBox platform](https://www.virtualbox.org/wiki/Downloads) should be installed on th host computer. This software runs under Linux, Solarix, Mac OSX and Windows operating systems.

2. A  sufficient amount of memory (RAM) to allocate two Gb to the RSAT-VM.

3. The VM initially occupies ~5Gb of hard disk space. By default, its hard drive is configured to grow dynamically when needed (for example when installing genomes), wih a maximal size of 18Gb.

------------------------------------

## Technical specifications

- *RSAT-VM version* : rsat-vm-2015-02
- *Linux distribution*: [ubuntu-14.04.1-desktop-amd64](http://releases.ubuntu.com/14.04.1/ubuntu-14.04.1-desktop-amd64.iso)

------------------------------------

## Downloading RSAT-VM

The Virtual Machine is available at [http://teaching.rsat.eu/](http://teaching.rsat.eu/), by clicking on the link "Download". 

------------------------------------

## Importing the VM

After having downloaded the appliance (rsa-vm-YYYY-MM.ova), open it with the Virtual Box application. 

Virtual Box proposes you to import the appliance, i.e. create an instance on which you will be able to work. 
Before importing the appliance, we recommend to check the option *Reinizialize MAC address of all network cards* at the bottom of the import dialog box. 

## Configuration


The RSAT-VM is provided in a ready-to-use mode.  However, you may need to adapt the configuration of your VirtualBox environment in order to obtain a correct behaviour of the VM.

###  VirtualBox host-only adapter

<div class="protocol">
1. Open the <span class="program">VirtualBox</span> program.

1. Open <span class="option">VirtualBox Preferences</span>.
<br>Click on the <span class="option">Network</span> option.
<br>Click on the tab <span class="option">Host-only Networks</span>. Check if a host-only adapter is already installed. If not, create a new one by clicking the <span class="option">**+**</span> icon on the right side (<a href="images/vb_prefs_host-only_network.png">snapshot</a>). <a href="images/vb_prefs_host-only_network.png"><img width="400" border="1" src="images/vb_prefs_host-only_network.png"></a>
	
1. Double-click on the host-only adapter to change its parameters.

1. In the <span class="option">Adaptor</span> tab, set the parameters as follows (<a href="images/vb_prefs_host-only_adaptor.png">snapshot</a>).	 
	    1. IPv4 Address: 192.168.56.1
	    1. IPv4 Network Mask: 255.255.255.0
	    1. IPv6 Address: (leave this field blank)
	    1. IPv5 Nework Mask length: 0
	   <a href="images/vb_prefs_host-only_network.png"><img width="400" border="1" src="images/vb_prefs_host-only_adaptor.png"></a>
	

1. In the tab <span class="option">DHCP Server</span>, set the parameters as follows (<a href="images/vb_prefs_host-only_dhcp.png">snapshot</a>).
	    1. Check the option <span class="option">Enable Server</span>
	    1. Server Address: 192.168.56.100
	    1. Server mask: 255.255.255.0
	    1. Lower Address Bound: 192.168.56.101
	    1. Upper Address Bound: 192.168.56.154
	 
<a href="images/vb_prefs_host-only_network.png"><img width="400" border="1" src="images/vb_prefs_host-only_dhcp.png"></a>

</div>


### Network settings for the guest machine


VirtualBox supports various ways to connect the guest (virtual machine) to the network.


#### Host-only network


This solution offers a good tradeoff between security and confort: your virtual machine (the guest) will be accessible only from your computer (the host).

<div class="protocol">
1. In the panel showing the available virtual machines, right-click on the RSAT-VM (rsat-vb-ub14d), open the <span class="option">Settings ...</span> dialog box.

1. In the tab <span class="option">Network</span>, select <span class="option">Adapter 1</span>, check <span class="option">Enable Network Adapter</span>, select <span class="option">Attached to: **Host-only Adapter**</span> (<a href="images/vm_settings_network_host-only.png">snapshot</a>).

1. In the pop-up menu besides the option <span class="option">Name</span>, select <span class="option">vboxnet0</span>.

<a href="images/vm_settings_network_host-only.png"><img width="400" border="1" src="images/vm_settings_network_host-only.png"></a>
</div>

#### NAT

**Note:** the host-only adapter will enable you to establish a connection (Web browsing, ssh connection) from the hosting operating system (the usual environment of your computer) to the guest system (the virtual machine). however, this adapter does not allow to connect the external world from the guest.

In parallel to the host-only adapter, we thus recommend to enable the second adapter and select <span class="option">NAT</span>.

#### Bridged network


Alternatively , for the sake of flexibility, you might consider to use a bridged network. The bridged adapter is the most convenient, because it establishes a bidirectional connection between your VM (the guest) and the network. Your guest RSAT Web server can thus be used from any other computer in your network. This configuration can typically be usd to make an RSAT server available for all people from the same lab or institute.

<div class="attention">
 **Attention!**  The bridged network makes your virtual
	machine visible for all the other computers of the local
	network the host machine (your PC). Check with your system
	administrator that this fits the local security
	requirements.
</div>

<div class="protocol">

1. In the panel showing the available virtual machines, right-click on the RSAT-VM (rsat-vb-ub14d), open the <span class="option">Settings ...</span> dialog box.
1. In the tab <span class="option">Network</span>, select <span class="option">Adapter 1</span>, check <span class="option">Enable Network Adapter</span>, select <span class="option">Attached to: Bridge   adapter</span>.
1. In the pop-up menu besides the option <span class="option">Name</span>, select an adapter depending on your local network configuration, e.g. <span class="option">Wi-fi (Airport)</span> (<a href="images/vm_settings_network_bridged.png">snapshot</a>).
	    <a href="images/vm_settings_network_bridged.png"><img width="400" border="1" src="images/vm_settings_network_bridged.png"></a>
	
</div>


------------------------------------

## Running RSAT-VM



### Starting the RSAT virtual machine

<div class="protocol">

1. In the left panel of <span class="program">VirtualBox</span>, select the  virual machine (rsat-vb-ub14d), and click on   the <span class="option">Start</span> icon.
 
<p class="tips">
At this stage, your RSAT VM should now be ready to be used from the Web interface. Assuming that you activated the host-only network as described abve, and that you only started one virtual machine), VirtualBox hould have assigned the first IP address of the range defined in the settings above: 192.168.56.101.</p>

</div>

### Using RSAT-VM as Web server

<div class="protocol">

1.  Open a connection to **<a target="_blank"
					href="http://192.168.56.101/">http://192.168.56.101/</a>**   in your web browser.
	 
<p class="tips">If the link does not work, it probably means that your network was not activated as described above. You will then need to obtain the IP address of your VM. Unfortunately, VirtualBox does not provide a direct way to know which IP address has been assigned to a VM. The only way we found to get this information is to <a href="#rsat-vm_login">log in in the VM</a>, open a terminal, and run the command <span class="command">/sbin/ifconfig</span>.</p>

</div>


### RSAT-VM log-in

You can log with the following parameters:

1. Username: vmuser
1. Password: tochng
 


<p class="tips">We intently chose an overly simple temporary password to ensure compatibility with AZERTY as well as QUERTY keyboards, but we recommend to use a safer password.</p>

<div class="attention">
 
1. At your first login, you will be prompted to change your password before anything else.

1. The user vmuser is sudoer. After login, you can thus become the master of your Virtual Machine, create new users, install packages, etc.
 
</div>

### Chosing the adequate keyboard for your computer

A small difficulty when distributing a VM is the large variety of keybords expected to be found on the users' computers. By default, we selected the standard British QWERTY keyboard.


On Ubuntu 14.04 server version, keybord configuration can be
  modified with following command.

<pre class="brush:bash;">
  sudo dpkg-reconfigure console-data
</pre>

For the desktop version, click on the <span class="option">Settings</span> icon <a href="images/icon_settings.png"><img width="40" border=1 src="images/icon_settings.png"></a>, then on the <a href="images/text-entry_options.png"><img width="50" src="images/icon_text-entry.png" border=1>, and check the keyboard.</a>.


### Connecting RSAT-VM in ssh

<pre class="brush:bash;">
	ssh vmuser@192.168.56.101
</pre>



