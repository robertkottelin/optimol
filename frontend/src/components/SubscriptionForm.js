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

const SubscriptionForm = ({ onSuccess, isMobile, isAuthenticated }) => {
  const stripe = useStripe();
  const elements = useElements();
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);
  const [error, setError] = useState("");
  const [paymentProcessingStep, setPaymentProcessingStep] = useState(null);
  const [requiresAction, setRequiresAction] = useState(false);
  const { currentUser, token, registerAndSubscribe, login } = useContext(AuthContext);
  
  // Form states for guest registration
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [confirmPassword, setConfirmPassword] = useState("");
  
  const apiBaseUrl = "/api";

  const handleSubmit = async (event) => {
    event.preventDefault();

    // Prevent multiple submissions
    if (isSubscribeLoading) {
      return;
    }

    // Reset states
    setError("");
    setPaymentProcessingStep("initializing");
    setIsSubscribeLoading(true);
    setRequiresAction(false);

    try {
      // Get card element
      const cardElement = elements.getElement(CardElement);
      if (!cardElement) {
        setError("Card information not found. Please refresh and try again.");
        setIsSubscribeLoading(false);
        return;
      }

      // Form validation
      if (!isAuthenticated) {
        if (!email || !password || !confirmPassword) {
          setError("All fields are required");
          setIsSubscribeLoading(false);
          return;
        }

        if (!/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(email)) {
          setError("Please enter a valid email address");
          setIsSubscribeLoading(false);
          return;
        }

        if (password !== confirmPassword) {
          setError("Passwords do not match");
          setIsSubscribeLoading(false);
          return;
        }

        if (password.length < 8) {
          setError("Password must be at least 8 characters long");
          setIsSubscribeLoading(false);
          return;
        }
      }

      setPaymentProcessingStep("creating_payment_method");
      
      // Create payment method with Stripe
      const paymentMethodResult = await stripe.createPaymentMethod({
        type: "card",
        card: cardElement,
        billing_details: {
          email: isAuthenticated ? currentUser.email : email,
          ...(postalCode && { address: { postal_code: postalCode } })
        }
      });

      if (paymentMethodResult.error) {
        // Handle specific card errors with actionable messages
        setError(`Payment method error: ${paymentMethodResult.error.message}`);
        setIsSubscribeLoading(false);
        return;
      }

      setPaymentProcessingStep("processing_payment");

      // Different handling based on authentication status
      if (isAuthenticated) {
        // User is already authenticated, use regular subscription flow
        const response = await axios.post(
          `${apiBaseUrl}/subscribe`, 
          { paymentMethodId: paymentMethodResult.paymentMethod.id },
          { 
            headers: {
              'Content-Type': 'application/json',
              'Authorization': `Bearer ${token}`
            }
          }
        );

        // Handle different response scenarios
        if (response.data.requires_action) {
          // Payment requires additional authentication (3D Secure)
          setPaymentProcessingStep("requires_authentication");
          setRequiresAction(true);

          const { error: confirmationError } = await stripe.confirmCardPayment(
            response.data.payment_intent_client_secret
          );

          if (confirmationError) {
            setError(`Payment authentication failed: ${confirmationError.message}`);
            setIsSubscribeLoading(false);
            return;
          }

          // Authentication succeeded
          setPaymentProcessingStep("payment_confirmed");
          onSuccess();
        } 
        else if (response.data.status === "pending") {
          // Payment is pending, but no action required
          setPaymentProcessingStep("payment_pending");
          
          // Confirm card payment to complete the transaction
          const { error: confirmError } = await stripe.confirmCardPayment(
            response.data.clientSecret
          );

          if (confirmError) {
            setError(`Payment confirmation error: ${confirmError.message}`);
            setIsSubscribeLoading(false);
            return;
          }

          setPaymentProcessingStep("payment_confirmed");
          onSuccess();
        }
        else if (response.data.success) {
          // Payment immediately succeeded
          setPaymentProcessingStep("subscription_active");
          onSuccess();
        } 
        else {
          // Unexpected response
          setError(`Subscription error: ${response.data.error || "Unknown error"}`);
          setIsSubscribeLoading(false);
        }
      } else {
        // Try logging in first to check if user already exists
        try {
          setPaymentProcessingStep("checking_existing_account");
          const loginResult = await login(email, password);
          
          if (loginResult.success) {
            // User already exists and credentials are correct
            // Send another request for subscription with the new token
            setPaymentProcessingStep("processing_existing_user");
            
            const subscribeResponse = await axios.post(
              `${apiBaseUrl}/subscribe`, 
              { paymentMethodId: paymentMethodResult.paymentMethod.id },
              { 
                headers: {
                  'Content-Type': 'application/json',
                  'Authorization': `Bearer ${loginResult.token}`
                }
              }
            );

            // Handle subscription response similar to authenticated branch
            if (subscribeResponse.data.requires_action) {
              setPaymentProcessingStep("requires_authentication");
              setRequiresAction(true);

              const { error: confirmationError } = await stripe.confirmCardPayment(
                subscribeResponse.data.payment_intent_client_secret
              );

              if (confirmationError) {
                setError(`Payment authentication failed: ${confirmationError.message}`);
                setIsSubscribeLoading(false);
                return;
              }

              setPaymentProcessingStep("payment_confirmed");
              onSuccess();
            } 
            else if (subscribeResponse.data.success) {
              setPaymentProcessingStep("subscription_active");
              onSuccess();
            } 
            else {
              setError(`Subscription error: ${subscribeResponse.data.error || "Unknown error"}`);
              setIsSubscribeLoading(false);
            }
            
            return;
          } 
          
          // Login failed, but might be because user doesn't exist
          // Continue with registration process
        } catch (loginErr) {
          // Silent handling - continue with registration flow
        }

        // Proceed with registration and subscription for new user
        setPaymentProcessingStep("creating_new_account");
        const result = await registerAndSubscribe(email, password, paymentMethodResult.paymentMethod.id);
        
        if (result.success) {
          // Check if we need to confirm the payment (SCA/3DS)
          if (result.clientSecret) {
            setPaymentProcessingStep("requires_authentication");
            
            const { error: confirmError } = await stripe.confirmCardPayment(result.clientSecret);

            if (confirmError) {
              setError(`Payment authentication failed: ${confirmError.message}`);
              setIsSubscribeLoading(false);
              return;
            }
          }

          setPaymentProcessingStep("subscription_active");
          onSuccess();  // Notify parent component
        } else {
          const errorMsg = result.error || "Registration and subscription failed";
          
          // Check for specific error types
          if (errorMsg.includes("email already registered")) {
            setError("An account with this email already exists. Please log in instead.");
          } else if (errorMsg.includes("declined")) {
            setError("Your card was declined. Please try a different payment method.");
          } else {
            setError(errorMsg);
          }
          
          setIsSubscribeLoading(false);
        }
      }
    } catch (err) {
      setPaymentProcessingStep(null);
      console.error("Error in subscription process:", err);
      
      // Extract error message from response if available
      let errorMessage = "Payment processing failed";
      
      if (err.response?.data?.error) {
        errorMessage = err.response.data.error;
      } else if (err.message) {
        errorMessage = err.message;
      }
      
      // Make error messages more user-friendly
      if (errorMessage.includes("rate limit") || errorMessage.includes("429")) {
        errorMessage = "Too many payment attempts. Please try again in a few minutes.";
      } else if (errorMessage.includes("network") || errorMessage.includes("connection")) {
        errorMessage = "Network error. Please check your internet connection and try again.";
      }
      
      setError(`Subscription error: ${errorMessage}`);
      setIsSubscribeLoading(false);
    }
  };

  // Helper for getting the payment step message
  const getPaymentStepMessage = () => {
    switch (paymentProcessingStep) {
      case "initializing": return "Initializing payment...";
      case "creating_payment_method": return "Processing your card...";
      case "processing_payment": return "Processing payment...";
      case "checking_existing_account": return "Checking account...";
      case "processing_existing_user": return "Preparing subscription...";
      case "creating_new_account": return "Creating your account...";
      case "requires_authentication": return "Verifying payment...";
      case "payment_pending": return "Confirming payment...";
      case "payment_confirmed": return "Payment confirmed, finalizing...";
      case "subscription_active": return "Subscription activated!";
      default: return "Processing...";
    }
  };

  return (
    <div className={`subscription-form glass ${isMobile ? 'mobile-smaller-padding' : ''}`}>
      <h3 className={`form-title ${isMobile ? 'mobile-stack' : ''}`}>
        <span className="premium-badge"><SparkleIcon /> PRO</span>
        <span>Unlock Advanced Molecular Optimization (Contact for Enterprise Subscription)</span>
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
        {isAuthenticated ? (
          // For authenticated users, just show their email
          <div className="input-group">
            <p className={isMobile ? 'mobile-smaller-text' : ''} style={{ marginBottom: '8px' }}>
              Subscribing as: <strong>{currentUser?.email}</strong>
            </p>
          </div>
        ) : (
          // For guest users, show registration form
          <>
            <div className="input-group">
              <label 
                htmlFor="email" 
                style={{ 
                  display: "block", 
                  marginBottom: "5px", 
                  fontSize: "14px", 
                  fontWeight: "500", 
                  color: "#94a3b8" 
                }}
              >
                Email Address
              </label>
              <input
                id="email"
                type="email"
                value={email}
                onChange={(e) => setEmail(e.target.value)}
                className="email-input"
                placeholder="your@email.com"
                required
                disabled={isSubscribeLoading}
              />
            </div>

            <div className="input-group">
              <label 
                htmlFor="password" 
                style={{ 
                  display: "block", 
                  marginBottom: "5px", 
                  fontSize: "14px", 
                  fontWeight: "500", 
                  color: "#94a3b8" 
                }}
              >
                Password
              </label>
              <input
                id="password"
                type="password"
                value={password}
                onChange={(e) => setPassword(e.target.value)}
                style={{ 
                  width: "100%", 
                  padding: "12px 12px 12px 12px", 
                  backgroundColor: "rgba(15, 23, 42, 0.7)", 
                  border: "1px solid rgba(255, 255, 255, 0.1)", 
                  borderRadius: "6px", 
                  color: "#f0f4f8", 
                  fontSize: "14px" 
                }}
                placeholder="••••••••"
                required
                disabled={isSubscribeLoading}
                minLength={8}
              />
            </div>

            <div className="input-group">
              <label 
                htmlFor="confirmPassword" 
                style={{ 
                  display: "block", 
                  marginBottom: "5px", 
                  fontSize: "14px", 
                  fontWeight: "500", 
                  color: "#94a3b8" 
                }}
              >
                Confirm Password
              </label>
              <input
                id="confirmPassword"
                type="password"
                value={confirmPassword}
                onChange={(e) => setConfirmPassword(e.target.value)}
                style={{ 
                  width: "100%", 
                  padding: "12px 12px 12px 12px", 
                  backgroundColor: "rgba(15, 23, 42, 0.7)", 
                  border: "1px solid rgba(255, 255, 255, 0.1)", 
                  borderRadius: "6px", 
                  color: "#f0f4f8", 
                  fontSize: "14px" 
                }}
                placeholder="••••••••"
                required
                disabled={isSubscribeLoading}
                minLength={8}
              />
            </div>
          </>
        )}
                
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
              hidePostalCode: true, // We're collecting it separately
            }}
            className="stripe-element"
            disabled={isSubscribeLoading}
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
        
        {requiresAction && (
          <div style={{ 
            color: "#f59e0b", 
            marginBottom: "10px", 
            padding: "8px", 
            borderRadius: "6px",
            backgroundColor: "rgba(245, 158, 11, 0.1)",
            border: "1px solid rgba(245, 158, 11, 0.2)",
            fontSize: isMobile ? "12px" : "14px"
          }}>
            Your card requires additional verification. A popup window may appear.
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
              {getPaymentStepMessage()}
            </>
          ) : (
            isAuthenticated ? "Subscribe Now" : "Register & Subscribe"
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