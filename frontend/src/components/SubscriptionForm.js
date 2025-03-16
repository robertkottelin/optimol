import React, { useState, useContext } from "react";
import { useStripe, useElements, CardElement } from "@stripe/react-stripe-js";
import axios from "axios";
import { AuthContext } from "../AuthContext";

// SVG icons for enhanced UI
const LockIcon = () => (
  <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <rect x="5" y="11" width="14" height="10" rx="2" stroke="currentColor" strokeWidth="2"/>
    <path d="M8 11V7C8 4.79086 9.79086 3 12 3C14.2091 3 16 4.79086 16 7V11" stroke="currentColor" strokeWidth="2"/>
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
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);
  const [error, setError] = useState("");
  const { currentUser, token } = useContext(AuthContext);
  
  const apiBaseUrl = "https://optimizemolecule.com";

  const handleSubmit = async (event) => {
    event.preventDefault();

    const cardElement = elements.getElement(CardElement);

    if (!cardElement) {
      setError("Card element not found");
      return;
    }

    setIsSubscribeLoading(true);
    setError("");

    try {
      // Create payment method with Stripe
      const { error: stripeError, paymentMethod } = await stripe.createPaymentMethod({
        type: "card",
        card: cardElement,
      });

      if (stripeError) {
        setError(`Payment method error: ${stripeError.message}`);
        setIsSubscribeLoading(false);
        return;
      }

      // Send to backend with JWT auth via Authorization header
      const response = await axios.post(
        `${apiBaseUrl}/subscribe`, 
        { paymentMethodId: paymentMethod.id },
        { 
          withCredentials: true,
          headers: {
            'Content-Type': 'application/json',
            'Authorization': token ? `Bearer ${token}` : ''
          }
        }
      );

      if (response.data.success) {
        const { clientSecret } = response.data;

        // Confirm subscription with Stripe
        const { error: confirmError } = await stripe.confirmCardPayment(clientSecret);

        if (confirmError) {
          setError(`Subscription error: ${confirmError.message}`);
        } else {
          console.log("Subscription successful");
          onSuccess();  // Notify parent component
        }
      } else {
        setError(`Subscription error: ${response.data.error || "Unknown error"}`);
      }
    } catch (err) {
      console.error("Error in subscription process:", err);
      setError(`Subscription error: ${err.response?.data?.error || err.message || "Unknown error"}`);
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
          <p className={isMobile ? 'mobile-smaller-text' : ''} style={{ marginBottom: '8px' }}>
            Subscribing as: <strong>{currentUser?.email}</strong>
          </p>
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
        
        {error && (
          <div style={{ 
            color: "#f43f5e", 
            marginBottom: "10px", 
            padding: "8px", 
            borderRadius: "6px",
            backgroundColor: "rgba(244, 63, 94, 0.1)",
            border: "1px solid rgba(244, 63, 94, 0.2)",
            fontSize: isMobile ? "12px" : "14px"
          }}>
            {error}
          </div>
        )}
        
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