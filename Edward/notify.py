# Script copied and modified from:
# https://www.cmu.edu/computing/services/comm-collab/email-calendar/exchange/how-to/imap-exchange.html
# https://towardsdatascience.com/notify-with-python-41b77d51657e

import smtplib
import socket
import os
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart

def message(subject="Python Notification", text="", img=None, attachment=None):
    # build message contents
    msg = MIMEMultipart()
    msg['Subject'] = subject  # add in the subject
    msg.attach(MIMEText(text))  # add text contents

    # check if we have anything given in the img parameter
    if img is not None:
        # if we do, we want to iterate through the images, so let's check that
        # what we have is actually a list
        if type(img) is not list:
            img = [img]  # if it isn't a list, make it one
        # now iterate through our list
        for one_img in img:
            img_data = open(one_img, 'rb').read()  # read the image binary data
            # attach the image data to MIMEMultipart using MIMEImage, we add
            # the given filename use os.basename
            msg.attach(MIMEImage(img_data, name=os.path.basename(one_img)))

    # we do the same for attachments as we did for images
    if attachment is not None:
        if type(attachment) is not list:
            attachment = [attachment]  # if it isn't a list, make it one
            
        for one_attachment in attachment:
            with open(one_attachment, 'rb') as f:
                # read in the attachment using MIMEApplication
                file = MIMEApplication(
                    f.read(),
                    name=os.path.basename(one_attachment)
                )
            # here we edit the attached file metadata
            file['Content-Disposition'] = f'attachment; filename="{os.path.basename(one_attachment)}"'
            msg.attach(file)  # finally, add the attachment to our message object
    return msg

# 'smtp.gmail.com'
# 'smtp.exchange.andrew.cmu.edu'

# If you want to use a gmail address instead of cmu, replace the server below. 
# PARAMETER

def send(msg, server='smtp.exchange.andrew.cmu.edu', port='587'):
    # contain following in try-except in case of momentary network errors
    try:
        # initialise connection to email server, the default is Gmail
        smtp = smtplib.SMTP(server, port)
        # this is the 'Extended Hello' command, essentially greeting our SMTP or ESMTP server
        smtp.ehlo()
        # this is the 'Start Transport Layer Security' command, tells the server we will
        # be communicating with TLS encryption
        smtp.starttls()
        
        # read email and password from file
        with open('emailANDpassword.txt', 'r') as fp:
            email, pwd = fp.read().split()
            
        # login to gmail server
        smtp.login(email, pwd)
        # send notification to self
        smtp.sendmail(email, email, msg.as_string())
        # disconnect from the server
        smtp.quit()
    except socket.gaierror:
        print("Network connection error, email not sent.")


# How to use:

# Create a file named emailANDpassword.txt formatted as follows:
# "egaremo@andrew.cmu.edu\npassword"

# notify.msg(
#     subject="Cashflow Model Completion",
#     text=(f'{len(model.output)} loans processed.\n'
#           f'Total runtime: {runtime}'),
#     img=[
#         '../vis/loan01_amortisation.png',
#         '../vis/loan07_amortisation.png',
#         '../vis/loan01_profit_and_loss.png',
#         '../vis/loan07_profit_and_loss.png'
#     ]
# )

# notify.send(msg)

if __name__ == '__main__':
    msg = message(
        subject="Notification script test",
        text=(f'Hello world!\n'
            f'How are you doing? Make sure to stretch those legs.')
    )
    send(msg)
    print("Email sent.")
