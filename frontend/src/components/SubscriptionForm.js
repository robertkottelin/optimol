import React, { useState } from "react";
import { useStripe, useElements, CardElement } from "@stripe/react-stripe-js";
import axios from "axios";

// SVG icons for enhanced UI
const LockIcon = () => (
  <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <rect x="5" y="11" width="14" height="10" rx="2" stroke="currentColor" strokeWidth="2"/>
    <path d="M8 11V7C8 4.79086 9.79086 3 12 3C14.2091 3 16 4.79086 16 7V11" stroke="currentColor" strokeWidth="2"/>
  </svg>
);

const EmailIcon = () => (
  <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <rect x="2" y="4" width="20" height="16" rx="2" stroke="currentColor" strokeWidth="2"/>
    <path d="M22 7L13.03 12.7C12.7213 12.8934 12.3643 13 12 13C11.6357 13 11.2787 12.8934 10.97 12.7L2 7" stroke="currentColor" strokeWidth="2"/>
  </svg>
);

const CheckIcon = () => (
  <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <path d="M5 12L10 17L19 8" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
  </svg>
);

const SparkleIcon = () => (
  <svg width="18" height="18" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <path d="M12 3L13.4323 8.86433H19.5H13.4323L12 3Z" fill="currentColor"/>
    <path d="M12 3L10.5677 8.86433H4.5H10.5677L12 3Z" fill="currentColor"/>
    <path d="M12 21L13.4323 15.1357H19.5H13.4323L12 21Z" fill="currentColor"/>
    <path d="M12 21L10.5677 15.1357H4.5H10.5677L12 21Z" fill="currentColor"/>
    <path d="M19.5 12L13.6357 13.4323V19.5V13.4323L19.5 12Z" fill="currentColor"/>
    <path d="M19.5 12L13.6357 10.5677V4.5V10.5677L19.5 12Z" fill="currentColor"/>
    <path d="M4.5 12L10.3643 13.4323V19.5V13.4323L4.5 12Z" fill="currentColor"/>
    <path d="M4.5 12L10.3643 10.5677V4.5V10.5677L4.5 12Z" fill="currentColor"/>
  </svg>
);

const SubscriptionForm = ({ onSuccess, isMobile }) => {
  const stripe = useStripe();
  const elements = useElements();
  const [email, setEmail] = useState("");
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);

  const handleSubmit = async (event) => {
    event.preventDefault();

    const cardElement = elements.getElement(CardElement);

    if (!email) {
      alert("Please enter your email.");
      return;
    }

    setIsSubscribeLoading(true);

    try {
      // Create payment method
      const { error, paymentMethod } = await stripe.createPaymentMethod({
        type: "card",
        card: cardElement,
        billing_details: { email },
      });

      if (error) {
        console.error("Payment method error:", error);
        alert(`Payment error: ${error.message}`);
        return;
      }

      // Send to backend
      const response = await axios.post("http://localhost:5000/subscribe", {
        email,
        paymentMethodId: paymentMethod.id,
      });

      if (response.data.success) {
        const { clientSecret } = response.data;

        // Confirm subscription
        const { error: confirmError } = await stripe.confirmCardPayment(
          clientSecret
        );

        if (!confirmError) {
          console.log("Subscription successful");
          onSuccess(email);
        } else {
          console.error("Error confirming subscription:", confirmError);
          alert(`Subscription error: ${confirmError.message}`);
        }
      } else {
        console.error("Subscription error:", response.data.error);
        alert(`Subscription error: ${response.data.error}`);
      }
    } catch (err) {
      console.error("Error in subscription process:", err);
      alert("An unexpected error occurred. Please try again.");
    } finally {
      setIsSubscribeLoading(false);
    }
  };

  return (
    <div className={`subscription-form glass ${isMobile ? 'mobile-smaller-padding' : ''}`}>
      <h3 className={`form-title ${isMobile ? 'mobile-stack' : ''}`}>
        <span className="premium-badge"><SparkleIcon /> PRO</span>
        <span>Unlock Advanced Molecular Optimization</span>
      </h3>
      
      <div className={`benefits-list ${isMobile ? 'mobile-smaller-padding' : ''}`}>
        <div className="benefit-item">
          <div className="check-icon"><CheckIcon /></div>
          <div className={isMobile ? 'mobile-smaller-text' : ''}>Up to 100,000 iterations for classical optimization</div>
        </div>
        <div className="benefit-item">
          <div className="check-icon"><CheckIcon /></div>
          <div className={isMobile ? 'mobile-smaller-text' : ''}>Up to 1,000 iterations for quantum optimization</div>
        </div>
        <div className="benefit-item">
          <div className="check-icon"><CheckIcon /></div>
          <div className={isMobile ? 'mobile-smaller-text' : ''}>Access to extended basis sets (6-311G, cc-pVDZ)</div>
        </div>
      </div>
      
      <form onSubmit={handleSubmit}>
        <div className="input-group">
          <div className="input-icon">
            <EmailIcon />
          </div>
          <input
            type="email"
            value={email}
            onChange={(e) => setEmail(e.target.value)}
            placeholder="Enter your email"
            required
            className="email-input"
          />
        </div>
        
        <div className="card-element-container">
          <div className="card-element-icon">
            <LockIcon />
          </div>
          <CardElement
            options={{
              style: {
                base: {
                  color: "#f0f4f8",
                  fontSize: isMobile ? "14px" : "16px",
                  fontFamily: "'Inter', sans-serif",
                  fontSmoothing: "antialiased",
                  "::placeholder": { color: "#94a3b8" },
                  iconColor: "#94a3b8",
                },
                invalid: { color: "#f43f5e" },
              },
              hidePostalCode: true,
            }}
            className="stripe-element"
          />
        </div>
        
        <div className={`secure-badge ${isMobile ? 'mobile-smaller-text' : ''}`}>
          <LockIcon /> Secure payment - $10/month
        </div>
        
        <button
          type="submit"
          disabled={!stripe || isSubscribeLoading}
          className={`subscribe-button ${isSubscribeLoading ? "loading" : ""}`}
        >
          {isSubscribeLoading ? (
            <>
              <span className="spinner"></span>
              Processing...
            </>
          ) : (
            "Subscribe Now"
          )}
        </button>
        
        <div className={`subscription-terms ${isMobile ? 'mobile-smaller-text' : ''}`}>
          Cancel anytime. Subscription renews monthly.
        </div>
      </form>
    </div>
  );
};

export default SubscriptionForm;