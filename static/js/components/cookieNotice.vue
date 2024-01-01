<template>
  <div class="cookie-popup">
    <div class="cookie-content">
      <p>
        RiboVision 2.0 uses essential cookies to enhance your experience on our website.
        By using our website, you consent to our use of cookies in accordance with our
        <a href="/privacy-policy">Privacy Policy</a>.
      </p>
      <button class="btn btn-primary" @click="acceptCookies">Accept</button>
    </div>
  </div>
</template>

<script>
export default {
  methods: {
    acceptCookies() {
      // Set a cookie or flag to indicate that the user has accepted cookies
      // You can use a library like js-cookie or localStorage for this purpose
      
      // After accepting cookies, emit the event to inform the parent component
      vm.cookiesAccepted = true;
    },
  },
};
</script>