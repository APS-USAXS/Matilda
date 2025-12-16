#!/usr/bin/env python3

'''
Docstring for CodeFragments.monitor_instrument
we will design code which will monitor USAXS instrument and send e-mail if it hangs

We need to watch few things:
1. epics pv usxLAX:state if USAXS is running.
2. if data fre being collected, need to figrue out how yet

'''



"""
Monitor an EPICS PV (usxLAX:state) and notify by e‑mail when it drops to 0.
"""

import time
import smtplib
from email.message import EmailMessage
import epics

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
PV_NAME        = "usxLAX:state"          # EPICS PV to monitor
CHECK_INTERVAL = 15                      # seconds
ALERT_EMAIL    = "ilavsky@anl.gov"       # recipient
SMTP_SERVER    = "smtp.gmail.com"        # change to your mail server
SMTP_PORT      = 587
SMTP_USER      = "your_email@gmail.com"  # change to your account
SMTP_PASS      = "your_app_password"     # change to your password / app‑specific password

# ------------------------------------------------------------------
# Helper: send e‑mail
# ------------------------------------------------------------------
def send_email(subject: str, body: str) -> None:
    """Send an e‑mail via the configured SMTP server."""
    msg = EmailMessage()
    msg["From"] = SMTP_USER
    msg["To"]   = ALERT_EMAIL
    msg["Subject"] = subject
    msg.set_content(body)

    with smtplib.SMTP(SMTP_SERVER, SMTP_PORT) as smtp:
        smtp.starttls()               # secure the connection
        smtp.login(SMTP_USER, SMTP_PASS)
        smtp.send_message(msg)

# ------------------------------------------------------------------
# Main monitoring loop
# ------------------------------------------------------------------
def main() -> None:
    pv = epics.PV(PV_NAME)

    if not pv.wait_for_connection(timeout=10):
        raise RuntimeError(f"Cannot connect to PV {PV_NAME}")

    print(f"Monitoring PV '{PV_NAME}'. Press Ctrl+C to stop.")

    try:
        while True:
            value = pv.get()          # read current PV value (blocks if no connection)
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {PV_NAME} = {value}")

            if value == 0:
                subject = f"⚠️ {PV_NAME} is 0"
                body = (f"The EPICS PV '{PV_NAME}' has value 0.\n\n"
                        "Please check the instrument as it may have stopped.")
                try:
                    send_email(subject, body)
                    print("E‑mail sent.")
                except Exception as e:
                    print(f"Failed to send e‑mail: {e}")

            time.sleep(CHECK_INTERVAL)

    except KeyboardInterrupt:
        print("\nStopped by user.")


if __name__ == "__main__":
    main()

